import scipy as sp
import numpy as np
from scipy.optimize import minimize, fmin, basinhopping, brute
import matplotlib.pyplot as plt
from itertools import cycle
from scipy.optimize import curve_fit

I = 2
skipGalinat = True
skipSimmons = False


def xAxis(x, aRe, aSt, aCa, C1, C2):
  Re = x[0]
  St = x[1]
  Ca = x[2]
  # works for I=3 and I=0
  #return Re ** aRe / St ** aSt / Ca ** aCa

  # works for I=2 
  return St ** aSt / Re ** aRe * Ca ** aCa

  # works for I=1
  #return Re ** aRe

def dependency(x, A, B, aRe, aSt, aCa, C1, C2):
  x = xAxis(x, aRe, aSt, aCa, C1, C2)
  return A * x + B

def getData():
    C = np.genfromtxt('validationData/coeffs_nondims.txt')
    if(skipSimmons):
      C = np.delete(C, 4, 0)
    C = C.T
    c = []
    if(skipGalinat):
      Re = C[4][4:]
      St = C[5][4:]
      Ca = C[6][4:]
      iMin = 4
      for i in range(4):
          c.append(C[i][4:])
    else:
      # Reynolds number in the orifice is twice as big as Re_pipe
      Re = C[4][1:4] * 2.0
      # St is 8 times as big
      St = C[5][1:4] * 8.0
      # Ca is 4 times as big
      Ca = C[6][1:4] * 4.0
      Re = np.append(Re, C[4][4:])
      St = np.append(St, C[5][4:])
      Ca = np.append(Ca, C[6][4:])
      for i in range(4):
          c.append(C[i][1:])


    return Re, St, Ca, np.array(c)

Re, St, Ca, C = getData()
#C[2] *= 1e14
#C[3] /= 1e12
X = [Re, St, Ca]
popt, pcov = curve_fit(dependency, X, C[I])
print('Error is ',  np.sqrt(np.diag(pcov)))
print('A: ', popt[0])
print('B: ', popt[1])
print('aRe: ', popt[2])
print('aSt: ', popt[3])
print('aCa: ', popt[4])
print('C1: ', popt[5])
print('C2: ', popt[6])

x = xAxis(X, *popt[2:])
y = dependency(X, *popt)
x_cont = np.linspace(min(x) / 2.0, max(x) * 2.0, 1000)

fig = plt.figure()
plt.plot(
    x, C[I], '+',
    linewidth=2, label="correlation")
plt.plot(
    x_cont, popt[0] * x_cont + popt[1],
    linewidth=2, label="correlation")
plt.xscale('log')
plt.yscale('log')
plt.show()
