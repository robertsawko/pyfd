from numpy import arange, zeros, sum, exp, trapz, piecewise, linspace
from scipy.integrate import odeint
from itertools import cycle
import matplotlib.pyplot as plt


class PBESolution:
    def pbe(
        self, N, t, xi, deltaXi,
        beta=lambda x, y: 2.0 / y, gamma=lambda x: x**2
    ):
        dim = N.shape[0]
        dNdt = - N * gamma(xi)
        for i in arange(dim):
            for j in arange(i + 1, dim):
                dNdt[i] += beta(xi[i], xi[j]) * gamma(xi[j]) * N[j] * deltaXi
        return dNdt

    def __init__(self, number_of_classes, t, xi0):
        N0 = zeros(number_of_classes)
        N0[-1] = 1
        self.xi = xi0 + xi0 * arange(number_of_classes)
        self.N = odeint(lambda NN, t: self.pbe(NN, t, self.xi, xi0), N0, t)


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


t = arange(0.0, 10.0, 0.001)
l = 1.0

pbe_solutions = dict(
    (n, PBESolution(n, t, l / n)) for n in [10, 20, 40, 80, 160]
)


totals = dict(
    (
        n,
        [sum(Ns) for Ns in pbe_solutions[n].N]
    ) for n in pbe_solutions
)

xi = linspace(0, l, 100)
Na = [zm_pure_breakup_total_number_solution(xi, time, l) for time in t]

fig = plt.figure()
ax = fig.gca()
linestyles = cycle(['-', '--', ':'])
for n in sorted(totals):
    ax.loglog(
        t, totals[n]/totals[n][0],
        linestyle=next(linestyles),
        label="MOC with N={0}".format(n))
ax.loglog(t, Na / Na[0], "-k", linewidth=2, label="Analytical")
ax.legend(loc='lower right', shadow=True)
ax.set_xlabel('t')
ax.set_ylabel('N/N0')
plt.show()

xi = linspace(0, l, 100, endpoint=False)
sol = zm_pure_breakup_pbe_solution(xi, t[-1], l)

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
for n in sorted(pbe_solutions):
    ax.semilogy(
        pbe_solutions[n].xi, pbe_solutions[n].N[-1] / (l/n), "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))
ax.semilogy(xi, sol, "-k", linewidth=2, label="Analytical")
ax.legend(loc='lower left', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')
plt.show()
