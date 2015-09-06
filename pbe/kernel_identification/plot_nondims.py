import scipy as sp
import numpy as np
from scipy.optimize import minimize, fmin, basinhopping, brute
import matplotlib.pyplot as plt
from itertools import cycle
import pylab

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
    C = data[3]
    #x = Re ** aRe * St ** aSt * Ca ** aCa
    if I == 2:
      #x = 0.062 * np.log(Re) * St ** 0.65776 * Ca ** 0.54667
      x = St ** 0.67935 / Re ** 0.9811 * Ca ** (-0.2008)
    elif I == 3:
      x = Re ** 1.6296 / St ** 0.5044 / Ca ** 0.2987
    else:
      x = Re

    ax.plot(
        x, C[I], '+', marker=next(markers),
        linewidth=2, label=names[name])

    xMin = min(xMin, np.amin(x))
    xMax = max(xMax, np.amax(x))
x = np.linspace(xMin, xMax)


if I == 2:
  #y = 1.7367e-14 + 1.6268e-12 * x
  y = - 3.1531e-14 + 3.33746e-09 * x
elif I == 1:
  y = 0.0355 + x * 0.0
elif I == 0:
  y = 0.00395 + x * 0.0
elif I == 3:
  y = 284255 * x


ax.grid()
if I > 1:
  ax.set_xscale('log')
  ax.set_yscale('log')

if I == 1 or I == 0:
  ax.set_xlabel(r'$Re$')
elif I == 2:
  ax.set_xlabel(r'$Re^{1.6296}St^{-0.5044} Ca^{-0.2987}$')
elif I == 3:
  ax.set_xlabel(r'$0.062 \log(Re) St^{0.65776} Ca^{0.54667}$')


ax.set_ylabel('C' + repr(I + 1))
ax.plot(
    x, y, color='black',
    linewidth=2)
plt.savefig('validationData/plots/C' + repr(I + 1) + '.pdf')

figLegend = pylab.figure(figsize = (4.5,1.5))
pylab.figlegend(*ax.get_legend_handles_labels(), loc = 'upper left')
plt.savefig('validationData/plots/legend.pdf')

plt.show()

