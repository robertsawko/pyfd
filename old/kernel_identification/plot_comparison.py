import numpy as np
from aux import plt, set_plt_params

set_plt_params()
fig = plt.figure()
ax = fig.gca()


dExp = np.loadtxt('validationData/comparison/dExp.txt') * 1e03
dNum = np.loadtxt('validationData/comparison/dNum.txt') * 1e03
x = np.linspace(0, 10)

ax.plot(dExp, dNum, 'ok')

line, = ax.plot(x, 1.1 * x, '--k', linewidth=0.5)
line, = ax.plot(x, 0.9 * x, '--k', linewidth=0.5)

line, = ax.plot(x, 1.2 * x, ':k', linewidth=0.25)
line, = ax.plot(x, 0.8 * x, ':k', linewidth=0.25)

line, = ax.plot(dExp, dExp, 'k')

ax.set_xlim(0, 1.5)
ax.set_ylim(0, 1.5)
ax.set_ylabel('d (numerical) [mm]')
ax.set_xlabel('d (experimental) [mm]')
fig.patch.set_alpha(0)
plt.savefig('validationData/plots/comparison.pgf', bbox_inches='tight')
