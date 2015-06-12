import scipy as sp
import numpy as np
from scipy.optimize import minimize, fmin, basinhopping, brute
import matplotlib.pyplot as plt
from itertools import cycle

aRe = -2
aSt = 1
aCa = 1
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
    c.append(C[0][1:4])
    c.append(C[1][1:4])
    c.append(C[2][1:4])
    c.append(C[3][1:4])

    return Re, St, Ca, c


def simmons():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    # Reynolds number in the orifice is twice as big as Re_pipe
    Re = C[4][4]
    # St is 8 times as big
    St = C[5][4]
    # Ca is 4 times as big
    Ca = C[6][4]
    c = []
    c.append(C[0][4])
    c.append(C[1][4])
    c.append(C[2][4])
    c.append(C[3][4])

    return Re, St, Ca, c


def coulaloglou():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    # Reynolds number in the orifice is twice as big as Re_pipe
    Re = C[4][5:19]
    # St is 8 times as big
    St = C[5][5:19]
    # Ca is 4 times as big
    Ca = C[6][5:19]
    c = []
    c.append(C[0][5:19])
    c.append(C[1][5:19])
    c.append(C[2][5:19])
    c.append(C[3][5:19])

    return Re, St, Ca, c


def angeli():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    # Reynolds number in the orifice is twice as big as Re_pipe
    Re = C[4][19:]
    # St is 8 times as big
    St = C[5][19:]
    # Ca is 4 times as big
    Ca = C[6][19:]
    c = []
    c.append(C[0][19:])
    c.append(C[1][19:])
    c.append(C[2][19:])
    c.append(C[3][19:])

    return Re, St, Ca, c

Re_sa, St_sa, Ca_sa, C_sa = simmons()
Re_g, St_g, Ca_g, C_g = galinat()
Re_c, St_c, Ca_c, C_c = coulaloglou()
Re_a, St_a, Ca_a, C_a = angeli()

fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])

x_g = Re_g ** aRe * St_g ** aSt * Ca_g ** aCa
x_c = Re_c ** aRe * St_c ** aSt * Ca_c ** aCa
x_sa = Re_sa ** aRe * St_sa ** aSt * Ca_sa ** aCa
x_a = Re_a ** aRe * St_a ** aSt * Ca_a ** aCa

if not skipGalinat:
    ax.plot(
        x_g, C_g[I], '+', marker=next(markers),
        linewidth=2, label="Galinat")

ax.plot(
    x_c, C_c[I], '+', marker=next(markers),
    linewidth=2, label="Coulaloglou")

ax.plot(
    x_sa, C_sa[I], '+', marker=next(markers),
    linewidth=2, label="Simmons")

ax.plot(
    x_a, C_a[I], '+', marker=next(markers),
    linewidth=2, label="Angeli")

xMin = min(np.amin(x_c), np.amin(x_a))
xMax = max(np.amax(x_c), np.amax(x_a))
if not skipGalinat:
    xMax = max(xMax, np.amax(x_g))
x = np.linspace(xMin, xMax)
# for C4 x = Re / St
#y = 3.0e08 * x
#for C3; x = St * Ca / Re ** 2
y = 1.0e-07 * x ** 0.5
# for C2; x = St
#y = np.exp(x / 14.0) / 40.0
#for C1; x = St
#y = np.exp(x / 14.0) / 400.0
ax.plot(
    x, y,
    linewidth=2, label="correlation")

ax.legend(loc='best', shadow=True)
ax.grid()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('St * Ca / Re^2')
ax.set_ylabel('C' + repr(I + 1))
plt.savefig('validationData/plots/C' + repr(I + 1) + '.pdf')
plt.show()
