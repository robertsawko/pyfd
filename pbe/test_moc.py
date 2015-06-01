from numpy import (
        exp, trapz, piecewise, arange, linspace, sqrt, zeros, sinh, array)
from numpy.testing import assert_almost_equal, assert_array_less
from moc import MOCSolution
from scipy.special import gamma


def check_error_convergence(L2_errors):
    # Testing convergence
    for k in arange(1, len(L2_errors)):
        assert_array_less(L2_errors[k], L2_errors[k - 1])

    # Testing convergence this will equal to less than 1% error
    assert_almost_equal(L2_errors[-1], 0.0, decimal=1)


def compare_with_analytical(
    time, grids, pbe_solutions, Na, pbe_analytical
):
    totals = dict(
        (
            n,
            array([sum(Ns) for Ns in pbe_solutions[n].N])
        ) for n in pbe_solutions
    )
    L2_total_errors = [
        L2_relative_error(time, totals[g], Na)
        for n, g in enumerate(grids)
    ]

    L2_pbe_errors = zeros(len(grids))
    for k, g in enumerate(grids):
        L2_pbe_errors[k] = L2_relative_error(
            pbe_solutions[g].xi,
            pbe_solutions[g].number_density()[-1],
            pbe_analytical(pbe_solutions[g].xi)
        )

    check_error_convergence(L2_total_errors)
    check_error_convergence(L2_pbe_errors)


def ziff_total_number_solution(x, t, l):
    """
    This is simply an integral of Ziff and McGrady
    """
    return exp(-t * l**2) \
        + trapz(2.0 * t * l * exp(-t * x**2), x=x)


def ziff_pbe_solution(x, t, l):
    """
    This is based on Equation 25 from Ziff and McGrady
    """
    return piecewise(
        x,
        [x < l, x == l, x > l],
        [
            lambda x: 2.0 * t * l * exp(-t * x**2),
            lambda x: exp(-t * x**2),
            lambda x: 0.0
        ]

    )


def scott_total_number_solution1(t):
    return exp(-t)


def scott_total_number_solution3(t, C=1, N0=1):
    T = C * N0 * t
    return 2.0 * N0 / (T + 2.0)


def scott_pbe_solution3(xi, t, C=1, N0=1, xi0=1):
    T = C * N0 * t
    x = xi / xi0
    #phi3 = 8.0 * exp(-2.0 * x) * sinh(2 * x * (T / (T + 2))**0.5) \
        #/ (T**0.5 * (T + 2)**1.5)
    phi3 = sum([
        (x * 2)**(2 * (k + 1)) / gamma(2 * (k + 1)) * (T/(T+2))**k
        for k in range(100)
        ]) * 4.0 * exp(- 2 * x) / (x * (T + 2)**2)
    return N0 / xi0 * phi3


def L2_relative_error(x, f, g):
    return sqrt(trapz((f - g)**2, x=x)) / sqrt(trapz(f**2, x=x))


def test_pure_binary_breakup():
    grids = [10, 20, 40, 80, 160]
    time = arange(0.0, 10.0, 0.001)
    l = 1.0

    pbe_solutions = dict()
    for g in grids:
        N0 = zeros(g)
        N0[-1] = 1
        pbe_solutions[g] = MOCSolution(
            N0, time, l / g,
            beta=lambda x, y: 2.0 / y,
            gamma=lambda x: x**2
        )

    # Calculate analytical solution to total numbers
    v = linspace(0, l, 100)
    Na = array([ziff_total_number_solution(v, t, l) for t in time])

    compare_with_analytical(
        time, grids, pbe_solutions, Na,
        lambda xi: ziff_pbe_solution(xi, time[-1], l)
    )


def test_pure_coalescence_constant():
    t = arange(0.0, 1, 0.01)
    vmax = 1e1
    v0 = 0.5
    N0 = 2
    grids = [10, 20, 40, 80, 160]
    C = 0.1

    pbe_solutions = dict()
    for g in grids:
        dv = vmax / g
        v = dv + dv * arange(g)
        Ninit = (N0 / v0) * (v / v0) * exp(-v / v0) * dv
        pbe_solutions[g] = MOCSolution(
            Ninit, t, dv,
            Q=lambda x, y: C
        )

    Na = scott_total_number_solution3(t, C=C, N0=N0)

    compare_with_analytical(
        t, grids, pbe_solutions, Na,
        lambda xi: scott_pbe_solution3(xi, t[-1], C=C, xi0=2 * v0, N0=N0)
    )
