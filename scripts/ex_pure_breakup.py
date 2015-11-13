from numpy import arange, linspace, array, piecewise
from itertools import cycle
from pyfd.pbe.moc import MOCSolution
from pyfd.tests.moc import ziff_total_number_solution, ziff_pbe_solution
import matplotlib.pyplot as plt

"""
Case setup based on:

    Ziff, R.M and McGrady, E.D
    "New solutions to the fragmentation equation", J. Phys A, vol. 18, 1985

    We are looking at case 4. from their paper with kernel F(x, y) = x + y.
    This corresponds to choosing:
        beta = 1.0/y,
        Gamma = y^2,
    for our kernels.
"""

grids = [10, 20, 40, 80, 160]
time = arange(0.0, 10.0, 0.001)
vmax = 1.0
pbe_solutions = dict()

for g in grids:
    threshold = vmax / g / 2

    # This is modelling Dirac's delta
    def N0init(x):
        return piecewise(
            x, [x < vmax - threshold, x >= vmax - threshold], [0, g / vmax])
    pbe_solutions[g] = MOCSolution(
        g, time, vmax / g, N0=N0init,
        beta=lambda x, y: 1.0 / y,
        gamma=lambda x: x**2
    )

totals = dict((n, pbe_solutions[n].total_numbers) for n in pbe_solutions)

v = linspace(0, vmax, 100)
Na = array([ziff_total_number_solution(v, t, vmax) for t in time])

fig = plt.figure()
ax = fig.gca()
linestyles = cycle(['-', '--', ':'])
for n in sorted(totals):
    ax.loglog(
        time, totals[n] / totals[n][0],
        linestyle=next(linestyles),
        label="MOC with N={0}".format(n))
ax.loglog(time, Na, "-k", linewidth=2, label="Analytical")
ax.legend(loc='lower right', shadow=True)
ax.set_xlabel('t')
ax.set_ylabel('N/N0')
plt.show()


fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
v = linspace(0, vmax, 10000)

for n in sorted(pbe_solutions):
    ax.loglog(
        pbe_solutions[n].xi, pbe_solutions[n].number_density[-1], "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))
ax.loglog(
    v, ziff_pbe_solution(v, time[-1], vmax), "-k",
    linewidth=2, label="Analytical $t=\infty$")
ax.legend(loc='upper right', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')
plt.show()
