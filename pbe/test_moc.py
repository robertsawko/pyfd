from numpy import exp, trapz, piecewise, arange, linspace, sqrt, zeros, sinh
from numpy.testing import assert_almost_equal, assert_array_less
from moc import MOCSolution
from scipy.special import gamma


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

    totals = dict(
        (
            n,
            [sum(Ns) for Ns in pbe_solutions[n].N]
        ) for n in pbe_solutions
    )
    v = linspace(0, l, 100)
    Na = [ziff_total_number_solution(v, t, l) for t in time]
    L2_total_errors = [
        L2_relative_error(time, totals[g]/totals[g][0], Na/Na[0])
        for n, g in enumerate(grids)
    ]

    L2_pbe_errors = zeros(len(grids))
    for k, g in enumerate(grids):
        L2_pbe_errors[k] = L2_relative_error(
            pbe_solutions[g].xi,
            pbe_solutions[g].N[-1] / (l / g),
            ziff_pbe_solution(pbe_solutions[g].xi, time[-1], l)
        )

    # Testing convergence
    for k in arange(1, len(grids)):
        assert_array_less(L2_total_errors[k], L2_total_errors[k - 1])
        assert_array_less(L2_pbe_errors[k], L2_pbe_errors[k - 1])

    # Testing convergence this will equal to less than 1% error
    assert_almost_equal(L2_total_errors[-1], 0.0, decimal=1)
    assert_almost_equal(L2_pbe_errors[-1], 0.0, decimal=1)
