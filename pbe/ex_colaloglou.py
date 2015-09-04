from numpy import arange, linspace, array, sqrt, exp, pi, zeros
from scipy.integrate import trapz
from itertools import cycle
from moc import MOCSolution
import matplotlib.pyplot as plt


"""
Case setup based on:

    Colaloglou and Tavlarides (1977)
    "", J. Chem Eng

"""

m_to_cm = 100
m3_to_cm3 = m_to_cm**3
kg_to_g = 1000
N_to_dyne = 10**5

D = 10  # [cm] Impeller diameter
Nstar = 4.16  # rps
# Water
muc = 1.0e-3 * kg_to_g / m_to_cm  # [g * cm^-1 s^-1]
rhoc = 1000.0 / m3_to_cm3 * kg_to_g
# Kerosene-dicholorebenzene
rhod = 972.0 / m3_to_cm3 * kg_to_g
sigma = 42.82e-03 / m_to_cm * N_to_dyne

phi = 0.15  # Holdup

M = 80
grids = [M]
time = arange(0.0, 100.0, 1)

mm3_to_cm3 = 0.1**3
vmax = 0.1 * mm3_to_cm3

# Feed distribution
v0 = 0.03 * mm3_to_cm3
sigma0 = 0.01 * mm3_to_cm3

# Feed
theta = 600
n0 = 0.01


def g(v, C1=0.336, C2=0.106):
    return C1 * v**(-2. / 9) * D**(2. / 3) * Nstar * \
        exp(- C2 * sigma / (rhod * v**(5. / 9) * D**(4. / 3) * Nstar**2))


def Q(v1, v2, C3=2.32e-6, C4=1.2e9):
    d_ratio = (v1**(1. / 3) * v**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))

    return C3 * \
        (v1**(2. / 3) + v2**(2. / 3)) * \
        (v1**(2. / 9) + v2**(2. / 9)) * \
        D**(2. / 3) * \
        Nstar * \
        exp(-C4 * muc * rhoc * D**2 / sigma**2 * Nstar**3 * d_ratio**4)


def beta(v1, v2):
    return 2. * 2.4 / v2 * exp(-4.5 * (2 * v1 - v2)**2 / (v2**2))

pbe_solutions = dict()
for m in grids:
    Ninit = zeros(m)
    Ninit[-1] = 1
    dv = vmax / m
    v = dv + dv * arange(m)
    A0 = phi / (v0 * sigma0 * sqrt(2 * pi)) * \
        exp(-(v - v0)**2 / (2 * sigma0**2)) * dv
    pbe_solutions[m] = MOCSolution(
        Ninit, time, vmax / m,
        beta=beta, gamma=g, Q=Q,
        theta=theta, n0=n0, A0=A0
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
    ax.plot(
        pbe_solutions[n].xi,
        pbe_solutions[n].number_density()[-1], "+-",
        # trapz(pbe_solutions[n].number_density()[-1], x=pbe_solutions[n].xi),
        label="MOC with N={0}".format(n))

ax.plot(
    pbe_solutions[M].xi,
    pbe_solutions[M].number_density()[0],  # /
    # trapz(pbe_solutions[M].number_density()[0], x=pbe_solutions[M].xi),
    label="Initial condition".format(n))

ax.legend(loc='upper right', shadow=True)
ax.set_xlabel('Particle volume')
ax.set_ylabel('Number density function')
plt.show()
