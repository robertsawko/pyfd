from numpy import arange, linspace, array, sqrt, exp, pi
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt


"""
Case setup based on:

    Colaloglou and Tavlarides (1977)
    "New solutions to the fragmentation equation", J. Chem Eng 

    We are looking at case 4. from their paper with kernel F(x, y) = x + y.
    This corresponds to choosing:
        beta = 2.0/y,
        Gamma = y^2,
    for our kernels.
"""

cm = 0.1
cm3 = cm**2
D = 0.1  # Impeller diameter
Nstar = 310.  # rpm
# Water
muc = 1.0e-3 * 1000 / cm
rhoc = 1000.0 / cm3 * 1000
# Kerosene-dicholorebenzene
rhod = 972.0 / cm3 * 1000
sigma = 42.82e-03 / cm

phi = 0.1  # Holdup

grids = [20, 40]
time = arange(0.0, 100.0, 0.001)
vmax = 1 * cm3

# Feed distribution
N0 = 100.0
v0 = 0.2
sigma0 = 0.01

# Feed
theta = 600
u0 = 0.05


def g(v, C1=0.4, C2=0.08):
    return C1 * v**(-2. / 9) * D**(2. / 3) * Nstar / (1 + phi) * \
        exp(- C2 * sigma * (1 + phi)**2 /
            (rhod * v**(5. / 9) * D**(4. / 3) * Nstar**2))


def Q(v1, v2, C3=2.8e-6, C4=1.83e9):
    d_ratio = (v1**(1. / 3) * v2**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))

    return C3 * \
        (v1**(2. / 3) + v2**(2. / 3)) * \
        (v1**(2. / 9) + v2**(2. / 9)) * \
        D**(2. / 3) * \
        Nstar / (1 + phi) * \
        exp(-C4 * muc * rhoc * D**2 / sigma**2 * Nstar**3 / (1 + phi)**3 *
            d_ratio**4)


def beta(v1, v2):
    return 2. * 2.4 / v1 * exp(-4.5 * (2 * v2 - v1)**2 / (v1**2))

pbe_solutions = dict()
for m in grids:
    dv = vmax / m
    v = dv + dv * arange(m)
    Ninit = N0 * phi / (v0 * sigma0 * sqrt(2 * pi)) *\
        exp(-(v - v0)**2 / (2 * sigma0**2)) * dv
    pbe_solutions[m] = MOCSolution(
        Ninit, time, vmax / m,
        beta=beta, gamma=g, Q=Q,
        theta=theta, u0=u0
    )

totals = dict(
    (
        n,
        array([sum(Ns) for Ns in pbe_solutions[n].N])
    ) for n in pbe_solutions
)

v = linspace(0, vmax, 100)

fig = plt.figure()
ax = fig.gca()
linestyles = cycle(['-', '--', '-.', ':'])
for n in sorted(totals):
    ax.loglog(
        time, totals[n] / totals[n][0],
        linestyle=next(linestyles),
        label='MOC with N={0}'.format(n))
ax.set_xlabel('Time $t$')
ax.set_ylabel('Total number of drops $N/N_0$')
# plt.show()


fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
v = linspace(0, vmax, 10000)

for n in sorted(pbe_solutions):
    ax.semilogy(
        pbe_solutions[n].xi, pbe_solutions[n].number_density()[-1],
        label="MOC with N={0}".format(n))

ax.semilogy(
    pbe_solutions[20].xi, pbe_solutions[20].number_density()[0],
    label="Initial guess".format(n))

ax.legend(loc='lower left', shadow=True)
ax.set_xlabel('Particle volume')
ax.set_ylabel('Number density function')
plt.show()
