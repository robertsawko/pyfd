from numpy import arange, linspace, zeros
from itertools import cycle
from moc import MOCSolution
from test_moc import blatz_and_tobolsky_pbe_solution
import matplotlib.pyplot as plt
from steadyStateSolver import SteadyStateSolution
from time import time

"""
Case setup based on:

    Blatz, B.J and Tobolsky, V.
    "Note on the Kinetic of systems manifesting simulaneous
    polymerization-depolymerization phenomena"

    B&T equations capture the evolution of polymer sized from a distribution of
    monomers. The kernels give constant coefficients for breakup and
    coalescence. Special care needs to be taken as the case is discrete. 

    NOTE: B&T equations are normalised with respect to initial number of
    monomers N0. Because of the difference with our formulation the coalescence
    rates have to be divided by N0 to keep the problems equivalent. This is
    because density function appears twice.
"""

t = arange(0.0, 10, 0.005)
t_ss = arange(0.0, 5.0, 1.0)
# Initial number of monomers
N0 = 10000
grids = [20, 60]
#grids = [20]
kc = 1.0
kb = 0.25
vmax = 100

pbe_solutions = dict()

start = time()
for g in grids:
    Ninit = zeros(g)
    Ninit[0] = N0
    pbe_solutions[g] = MOCSolution(
        Ninit, t, 1.0,
        # Dividing coalescence coefficient by the number of monomers to make
        # formulations equivalent
        Q=lambda x, y: kc / N0,
        # Guard has to be imlemented in order to avoid division by zero
        beta=lambda x, y: 2.0 / max([y - 1.0, 1e-6]),
        gamma=lambda x: kb * (x - 1.0)
    )
print "Time elepsed for full PBE solution: ", time() - start

Ss = dict()
start = time()
for g in grids:
    Ninit = zeros(g)
    Ninit[0] = N0
    Ss[g] = SteadyStateSolution(
        Ninit, t, 1.0,
        Q=lambda x, y: kc / N0,
        # Guard has to be imlemented in order to avoid division by zero
        beta=lambda x, y: 2.0 / max([y - 1.0, 1e-6]),
        gamma=lambda x: kb * (x - 1.0)
    )
    # set time used to advance the solution before checking the current error
    Ss[g].setTime(t_ss)
    # set minimum and maximum number of iterations
    Ss[g].setIter(1, 200)
    # how much the error must decrease
    Ss[g].setError(1e-03)
    Ss[g].setMassTol(0.7, 1.3)
    Ss[g].solve()
print "Time elepsed for steady state PBE solution: ", time() - start

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
v = linspace(0.9, vmax, 10000)

ax.set_ylim([1e-10, 10e3])

for n in sorted(Ss):
    ax.loglog(
        Ss[n].xi, Ss[n].solution, "+",
        marker=next(markers),
        label="Steady state solution ({0} classes)".format(n))

for n in sorted(pbe_solutions):
    ax.loglog(
        pbe_solutions[n].xi, pbe_solutions[n].N[-1], "+",
        marker=next(markers),
        label="MOC with ({0} classes)".format(n))
ax.loglog(
    v, N0 * blatz_and_tobolsky_pbe_solution(v, t[-1], kc, kb), "-k",
    linewidth=2, label="Analytical $t=\infty$")
ax.legend(loc='best', shadow=True)
ax.set_xlabel('Volume')
ax.set_ylabel('N/N0')
plt.savefig('comparison.pdf')
plt.show()
