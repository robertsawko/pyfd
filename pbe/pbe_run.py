from numpy import arange, zeros, sum, exp, linspace
from itertools import cycle
from moc import MOCSolution
from test_moc import scott_pure_coalescence_total_number_solution
import matplotlib.pyplot as plt


t = arange(0.0, 0.001, 0.0001)
l = 1e3

pbe_solutions = dict()
for g in [10, 20]:
    v0 = l / g
    v = v0 + v0 * arange(g)
    N0 = exp(-v)
    pbe_solutions[g] = MOCSolution(
        N0, t, l / g,
        Q=lambda x, y: x + y
    )

totals = dict(
    (
        n,
        [sum(Ns) for Ns in pbe_solutions[n].N]
    ) for n in pbe_solutions
)

Na = [scott_pure_coalescence_total_number_solution(v, time, l) for time in t]

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

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
for n in sorted(pbe_solutions):
    ax.semilogy(
        pbe_solutions[n].xi, pbe_solutions[n].N[-1] / (l/n), "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))
# ax.semilogy(xi, sol, "-k", linewidth=2, label="Analytical")
ax.legend(loc='lower left', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')
plt.show()
