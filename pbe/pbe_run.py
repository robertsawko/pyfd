from numpy import arange, sum, exp, linspace
from itertools import cycle
from moc import MOCSolution
from test_moc import scott_total_number_solution3
from test_moc import scott_pbe_solution3
import matplotlib.pyplot as plt


t = arange(0.0, 1, 0.01)
vmax = 1e1
v0 = 0.5
N0 = 2
grids = [40, 80, 160]
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

#totals = dict(
    #(
        #n,
        #[sum(Ns) for Ns in pbe_solutions[n].N]
    #) for n in pbe_solutions
#)

#Na = scott_total_number_solution3(t)

#fig = plt.figure()
#ax = fig.gca()
#linestyles = cycle(['-', '--', ':'])
#for n in sorted(totals):
    #ax.loglog(
        #t, totals[n]/totals[n][0],
        #linestyle=next(linestyles),
        #label="MOC with N={0}".format(n))
#ax.loglog(t, Na, "-k", linewidth=2, label="Analytical")
#ax.legend(loc='lower right', shadow=True)
#ax.set_xlabel('t')
#ax.set_ylabel('N/N0')
#plt.show()

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
#xi = pbe_solutions[grids[-1]].xi
v = linspace(0, vmax, 100000)

ax.set_xlim([0, 5])
ax.set_ylim([1e-10, 1.5])
for n in sorted(pbe_solutions):
    ax.plot(
        pbe_solutions[n].xi, pbe_solutions[n].N[-1] / (vmax / n), "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))
ax.plot(
    v, N0/v0 * v/v0 * exp(-v/v0), "--k",
    linewidth=2, label="IC")
ax.plot(
    v, scott_pbe_solution3(v, t[-1], C=C, xi0=2*v0, N0=N0), "-k",
    linewidth=2, label="Analytical")
ax.legend(loc='upper right', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')
plt.show()
