from numpy import arange, sqrt, pi
from itertools import cycle
from fluid import fluid
import numpy as np
from aux import set_plt_params, plt

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
one_over_theta = 0.005
for i in np.arange(Q.shape[0]):
    if i != (Q.shape[0] - 1):
        for j in arange(Q.shape[0]):
            Q[i] = N[j] * F.Q(v[i], v[j])
    escape[i] = one_over_theta

set_plt_params(relative_fig_width=0.7)
fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
linestyles = cycle(['-', '--', ':', '-.'])

ax = fig.gca()
ax.plot(v / v0, gamma, label="breakup")
ax.plot(v / v0, Q, label="coalescence")
ax.plot(v / v0, escape, linestyle='--', label="escape", color='black')


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
    x / v0, y, 's', color='black')
ax.text(x[0] / v0, y[0] * 1.15, '3')
ax.text(x[1] / v0, y[1] * 1.15, '1')
ax.text(x[2] / v0, y[2] * 1.15, '2')

ax.legend(loc='best')
ax.set_xlim(0.5, 1.7)
ax.set_ylim(0.0, 0.013)
ax.set_xlabel(r'$v/v_0$')
ax.set_ylabel('rate [1/s]')
plt.yticks([one_over_theta], [r"$\frac{1}{\theta}$"])
plt.tick_params(
    axis='x', which='both', bottom='off', top='off', labelbottom='off')
plt.savefig('validationData/plots/rates.pgf', bbox_inches='tight')
plt.show()
