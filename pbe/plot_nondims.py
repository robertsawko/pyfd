import scipy as sp
import numpy as np
from scipy.optimize import minimize, fmin, basinhopping, brute
import matplotlib.pyplot as plt
from itertools import cycle

aRe = 1
aSt = 2
aCa = 0
I = 2
skipGalinat = True


def galinat():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    # Reynolds number in the orifice is twice as big as Re_pipe
    Re = C[4][1:4] * 2.0
    # St is 8 times as big
    St = C[5][1:4] * 8.0
    # Ca is 4 times as big
    Ca = C[6][1:4] * 4.0
    c = []
    for i in range(4):
        c.append(C[i][1:4])

    return Re, St, Ca, c


def simmons():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][4]
    St = C[5][4]
    Ca = C[6][4]
    c = []
    for i in range(4):
        c.append(C[i][4])

    return Re, St, Ca, c


def coulaloglou():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][5:19]
    St = C[5][5:19]
    Ca = C[6][5:19]
    c = []
    for i in range(4):
        c.append(C[i][5:19])

    return Re, St, Ca, c


def angeli():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][19:24]
    St = C[5][19:24]
    Ca = C[6][19:24]
    c = []
    for i in range(4):
        c.append(C[i][19:24])

    return Re, St, Ca, c


def karabelas():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][24:]
    St = C[5][24:]
    Ca = C[6][24:]
    c = []
    for i in range(4):
      c.append(C[i][24:])

    return Re, St, Ca, c


nondims = dict()
nondims['simmons'] = (simmons())
nondims['galinat'] = (galinat())
nondims['coulaloglou'] = (coulaloglou())
nondims['angeli'] = (angeli())
nondims['karabelas'] = (karabelas())

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])

xMin = 1e200
xMax = -1
for name in nondims:
    data = nondims[name]
    Re = data[0]
    St = data[1]
    Ca = data[2]
    C = data[3]
    x = Re ** aRe * St ** aSt * Ca ** aCa
    #x = Re / St ** 2

    if skipGalinat and name == 'galinat':
        continue
    ax.plot(
        x, C[I], '+', marker=next(markers),
        linewidth=2, label=name)

    xMin = min(xMin, np.amin(x))
    xMax = max(xMax, np.amax(x))
x = np.linspace(xMin, xMax)

# for C4 x = Re / St
#y = 3.0e08 * x
#for C3; x = St * Ca / Re ** 2
#y = 1.0e-07 * x ** 0.5
# for C2; x = St
y = np.exp(x / 14.0) / 40.0
#for C1; x = St
#y = np.exp(x / 14.0) / 400.0
#ax.plot(
    #x, y,
    #linewidth=2, label="correlation")

ax.legend(loc='best', shadow=True)
ax.grid()
ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlabel('St * Ca / Re^2')
ax.set_ylabel('C' + repr(I + 1))
plt.savefig('validationData/plots/C' + repr(I + 1) + '.pdf')
plt.show()
