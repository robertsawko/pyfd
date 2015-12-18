'''
Basset Bousinesq Oseen Equation

For now a simple code for calculating drag+virtual+gravity
'''

from numpy import pi, linspace
from scipy.integrate import odeint

dp = 1
mu = 1
rhof = 1
g = -9.81
rhop = 1.1835  # This should give -0.1 velocity

transient_coeff = pi / 6.0 * rhop * dp**3
buoyancy = pi / 6.0 * (rhop - rhof) * dp**3 * g
drag_coeff = 3 * pi * mu * dp
added_mass_coeff = pi / 12 * rhof * dp**3

lhs_coeff = transient_coeff + added_mass_coeff


def RHS(Up, t):
    return (-drag_coeff * Up + buoyancy) / lhs_coeff

t = linspace(0, 1, 1000)
result = odeint(RHS, 0.0, t)


from matplotlib.pyplot import figure

fig = figure()
ax = fig.gca()
ax.plot(t, result)
fig.show()
