import scipy as sp
import numpy as np
from scipy.optimize import minimize, fmin, basinhopping, brute
import matplotlib.pyplot as plt
from itertools import cycle
import pylab

I = 3
skipGalinat = True


def galinat():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    # Reynolds number in the orifice is twice as big as Re_pipe
    Re = C[4][1:4] * 2.0
    # St is 8 times as big
    St = C[5][1:4] * 8.0
    # Ca is 4 times as big
    Ca = C[6][1:4] * 4.0
    We = C[7][1:4]
    c = []
    for i in range(4):
        c.append(C[i][1:4])

    return Re, St, Ca, We, c


def simmons():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][4]
    St = C[5][4]
    Ca = C[6][4]
    We = C[7][4]
    c = []
    for i in range(4):
        c.append(C[i][4])

    return Re, St, Ca, We, c


def coulaloglou():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][5:19]
    St = C[5][5:19]
    Ca = C[6][5:19]
    We = C[7][5:19]
    c = []
    for i in range(4):
        c.append(C[i][5:19])

    return Re, St, Ca, We, c


def angeli():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][19:24]
    St = C[5][19:24]
    Ca = C[6][19:24]
    We = C[7][19:24]
    c = []
    for i in range(4):
        c.append(C[i][19:24])

    return Re, St, Ca, We, c


def karabelas():
    C = np.genfromtxt('validationData/coeffs_nondims.txt').T
    Re = C[4][24:]
    St = C[5][24:]
    Ca = C[6][24:]
    We = C[7][24:]
    c = []
    for i in range(4):
        c.append(C[i][24:])

    return Re, St, Ca, We, c


nondims = dict()
names = dict()

nondims['simmons'] = (simmons())
names['simmons'] = ('Simmons and Azzopardi (2001)')
if not skipGalinat:
    nondims['galinat'] = (galinat())
    names.append('Galinat et al. (2005)')
nondims['coulaloglou'] = (coulaloglou())
names['coulaloglou'] = ('Coulaloglou and Tavlarides (1977)')
nondims['angeli'] = (angeli())
names['angeli'] = ('Angeli and Hewitt (2000)')
nondims['karabelas'] = (karabelas())
names['karabelas'] = ('Karabelas (1978)')

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
    We = data[3]
    C = data[4]
    #x = Re ** aRe * St ** aSt * Ca ** aCa
    if I == 2:
        aCa = 1.3927752887548224
        aSt = 0.71738786775149754
        c1 = 3.5541424951605571
        c2 = 0.63108186158198221
        x = St ** aSt * (Ca ** aCa * c1 + c2)
    elif I == 3:
        aWe = 0.18740863508646571
        aSt = -0.6257571490782956
        x = Re * St ** aSt / We ** aWe
    elif I == 0:
        aWe = 1.1466914631499581
        aCa = -0.14907380395843822
        x = We ** aWe / Ca ** aCa
    elif I == 1:
        aWe = 1.6170490368507624
        aRe = -0.2238711406747752
        x = We ** aWe * Re ** aRe

    ax.plot(
        x, C[I], '+', marker=next(markers),
        linewidth=2, label=names[name])

    xMin = min(xMin, np.amin(x))
    xMax = max(xMax, np.amax(x))

x = np.linspace(xMin, xMax)
if I == 0:
    A = 0.11371549951931224
    B = 0.0020814470588332158
elif I == 1:
    A = 31.563299222959419
    B = 0.019272087314426867
elif I == 2:
    A = 1.6428900960606153e-13
    B = 6.1228384797660674e-15
elif I == 3:
    A = 513765273.33555937
    B = 0.0
y = A * x + B


ax.grid()
if I > 1:
    ax.set_xscale('log')
    ax.set_yscale('log')

if I == 0:
    # x = We ** aWe / Ca ** aCa
    ax.set_xlabel(
        r'$We^{' + format(aWe, '.3f') + '}$'
        + r'$Ca^{' + format(-aCa, '.3f') + '}$'
        )
if I == 1:
    # x = We ** aWe * Re ** aRe
    ax.set_xlabel(
        r'$We^{' + format(aWe, '.3f') + '}$'
        + r'$Re^{' + format(aRe, '.3f') + '}$'
        )
elif I == 2:
    # x = St ** aSt * (Ca ** aCa * c1 + c2)
    ax.set_xlabel(
        r'$St^{' + format(aSt, '.3f') + '}$'
        + r'$(' + format(c1, '.3f') + '$'
        r'$Ca^{' + format(aCa, '.3f') + '}$'
        r'$+' + format(c2, '.3f') + ')}$'
        )
elif I == 3:
    # x = Re * St ** aSt / We ** aWe
    ax.set_xlabel(
        r'$St^{' + format(aSt, '.3f') + '}$'
        + r'$We^{' + format(-aWe, '.3f') + 'Re}$'
        )


ax.set_ylabel('C' + repr(I + 1))
ax.plot(
    x, y, color='black',
    linewidth=2)
plt.savefig('validationData/plots/C' + repr(I + 1) + '.pdf')

figLegend = pylab.figure(figsize=(4.5, 1.5))
pylab.figlegend(*ax.get_legend_handles_labels(), loc='upper left')
plt.savefig('validationData/plots/legend.pdf')

plt.show()