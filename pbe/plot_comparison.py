import scipy as sp
import numpy as np
from scipy.optimize import minimize, fmin, basinhopping, brute
import matplotlib.pyplot as plt
from itertools import cycle

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])

dExp = np.loadtxt('validationData/comparison/dExp.txt') * 1e03
dNum = np.loadtxt('validationData/comparison/dNum.txt') * 1e03

#ax.legend(loc='best', shadow=True)
ax.grid()

ax.plot(
    dExp, dNum, '+', marker=next(markers), color='black',
    linewidth=2)

line, = ax.plot(
    dExp, 1.1 * dExp, linewidth=1, color='black', linestyle=':')
line.set_dashes([6, 20])
line, = ax.plot(
    dExp, 0.9 * dExp, linewidth=1, color='black', linestyle=':')
line.set_dashes([6, 20])

line, = ax.plot(
    dExp, 1.2 * dExp, linewidth=1, color='black', linestyle='--')
line.set_dashes([1, 10])
line, = ax.plot(
    dExp, 0.8 * dExp, linewidth=1, color='black', linestyle='--')
line.set_dashes([1, 10])

line, = ax.plot(
    dExp, dExp, linewidth=1, color='black')

ax.set_xlim(0, 1.5)
ax.set_ylim(0, 1.5)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylabel('d (numerical) [mm]')
ax.set_xlabel('d (experimental) [mm]')
plt.savefig('validationData/plots/comparison.pdf')
plt.show()
