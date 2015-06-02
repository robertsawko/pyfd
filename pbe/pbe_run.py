from numpy import arange, linspace, sqrt, zeros
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt

time = arange(0.0, 10, 0.005)
vmax = 1e4
v0 = 0.5
N0 = 10000
grids = [20, 40, 80, 160]
kc = 1.0
kb = 0.26

pbe_solutions = dict()

for g in grids:
    Ninit = zeros(g)
    Ninit[0] = N0
    pbe_solutions[g] = MOCSolution(
        Ninit, time, 1.0,
        Q=lambda x, y: kc / N0,
        beta=lambda x, y: 2.0 / max([y - 1.0, 1e-6]),
        gamma=lambda x: kb * (x - 1.0)
    )

from numpy import tanh
K = 1.0 + kb / kc
pinf = 1.0 / (K + sqrt(K**2 - 1))
p = 1.0 / (K + sqrt(K**2 - 1) * 1.0 / tanh(kc*time[-1] / 2.0 * sqrt(K**2 - 1)))


fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
v = linspace(0.9, vmax, 10000)

ax.set_ylim([1e-10, 10e4])
for n in sorted(pbe_solutions):
    ax.loglog(
        pbe_solutions[n].xi, pbe_solutions[n].N[-1], "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))
ax.loglog(
    v,  N0 * p**(v - 1.0) * (1.0 - p)**2, "--r",
    linewidth=2, label="Analytica $t={0}$".format(time[-1]))
ax.loglog(
    v, N0 * pinf**(v - 1.0) * (1.0 - pinf)**2, "-k",
    linewidth=2, label="Analytical $t=\infty$")
ax.legend(loc='upper right', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')
plt.show()
