from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from steadyStateSolver import SteadyStateSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np
from scipy.optimize import minimize, fmin, fmin_powell
from sys import argv

F = fluid('coulaloglou', 0, 0)
F.C = np.array([0.00299, 0.0256, 6.27e-11, 1.2e14])

d0 = 3e-04
v0 = np.pi * d0 ** 3 / 6.0
s0 = v0 / 2.0
v = np.arange(v0 / 10.0, 3.0 * v0, 0.01 * v0)
N = F.alpha / v0 * F.V\
    * 1.0 / s0 / sqrt(2.0 * pi)\
    * np.exp(- (v - v0) ** 2 / 2 / s0 ** 2)

gamma = F.gamma(v)
Q = np.zeros(v.shape)
escape = np.zeros(v.shape)

for i in np.arange(Q.shape[0]):
    if i != (Q.shape[0] - 1):
        for j in arange(Q.shape[0]):
            Q[i] = N[j] * F.Q(v[i], v[j])
    escape[i] = 0.005

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
linestyles = cycle(['-', '--', ':', '-.'])

ax = fig.gca()
ax.plot(
    v / v0, gamma, linewidth=2,
    label="breakup rate")
ax.plot(
    v / v0, Q, linewidth=2,
    label="coalescence rate")
ax.plot(
    v / v0, escape, linewidth=2, linestyle='--',
    label="escape frequency", color='black')


x = np.zeros(3)
y = np.zeros(3)
coalBreak_val = min(abs(Q - gamma)[:-1])
coalBreak = np.where(abs(Q - gamma) == coalBreak_val)
x[0] = v[coalBreak]
y[0] = Q[coalBreak]

escapeBreak_val = min(abs(escape - gamma)[:-1])
escapeBreak = np.where(abs(escape - gamma) == escapeBreak_val)
x[1] = v[escapeBreak]
y[1] = gamma[escapeBreak]

escapeCoal_val = min(abs(escape - Q)[:-1])
escapeCoal = np.where(abs(escape - Q) == escapeCoal_val)
x[2] = v[escapeCoal]
y[2] = Q[escapeCoal]

ax.plot(
    x / v0, y, 's', markersize=10, color='black')
ax.text(x[0] / v0, y[0] * 1.15, '3', fontsize=16)
ax.text(x[1] / v0, y[1] * 1.15, '1', fontsize=16)
ax.text(x[2] / v0, y[2] * 1.15, '2', fontsize=16)

ax.legend(loc='best')
ax.set_xlim(0.5, 1.7)
ax.set_ylim(0.0, 0.013)
ax.set_xlabel(r'$v/v_0$')
ax.set_ylabel('rate [1/s]')
plt.savefig('validationData/plots/rates.pdf')
plt.show()
