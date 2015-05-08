from numpy import exp, trapz, piecewise, arange, linspace, sqrt, zeros
from numpy.testing import assert_almost_equal, assert_array_less
from moc import MOCSolution


def zm_pure_breakup_total_number_solution(x, t, l):
    """
    This is simply an integral of Ziff and McGrady
    """
    return exp(-t * l**2) \
        + trapz(2.0 * t * l * exp(-t * x**2), x=x)


def zm_pure_breakup_pbe_solution(x, t, l):
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


def L2_relative_error(x, f, g):
    return sqrt(trapz((f - g)**2, x=x)) / sqrt(trapz(f**2, x=x))


def test_pure_binary_breakup():
    grids = [10, 20, 40, 80, 160]
    time = arange(0.0, 10.0, 0.001)
    l = 1.0

    pbe_solutions = dict(
        (n, MOCSolution(n, time, l / n)) for n in grids
    )

    totals = dict(
        (
            n,
            [sum(Ns) for Ns in pbe_solutions[n].N]
        ) for n in pbe_solutions
    )
    v = linspace(0, l, 100)
    Na = [zm_pure_breakup_total_number_solution(v, t, l) for t in time]
    L2_total_errors = [
        L2_relative_error(time, totals[g]/totals[g][0], Na/Na[0])
        for n, g in enumerate(grids)
    ]

    L2_pbe_errors = zeros(len(grids))
    for k, g in enumerate(grids):
        L2_pbe_errors[k] = L2_relative_error(
            pbe_solutions[g].xi,
            pbe_solutions[g].N[-1] / (l / g),
            zm_pure_breakup_pbe_solution(pbe_solutions[g].xi, time[-1], l)
        )

    # Testing convergence
    for k in arange(1, len(grids)):
        assert_array_less(L2_total_errors[k], L2_total_errors[k - 1])
        assert_array_less(L2_pbe_errors[k], L2_pbe_errors[k - 1])

    # Testing convergence this will equal to less than 1% error
    assert_almost_equal(L2_total_errors[-1], 0.0, decimal=1)
    assert_almost_equal(L2_pbe_errors[-1], 0.0, decimal=1)
