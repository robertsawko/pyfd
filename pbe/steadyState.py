from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from steadyStateSolver import SteadyStateSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np
from scipy.optimize import minimize, fmin, fmin_powell, differential_evolution
from sys import argv


def case_error(C):
    pbe_solutions = dict()
    F = fluid(argv[1], (int)(argv[2]), (int)(argv[3]))

    F.C1 = C[0]
    F.C2 = C[1]
    F.C3 = C[2]
    F.C4 = C[3]
    print '--------------'
    print 'constants: ', F.C1, F.C2, F.C3, F.C4

    Er = 0.0
    for i in np.arange(2):
      t = F.timeRange
      v0 = F.v0[i]
      s0 = F.s0[i]
      vmax = F.vMax[i]
      g = F.numberOfClasses
      print "initial d and v: ", i, F.d0[i], v0

      dv = vmax / g
      v = dv + dv * arange(g)
      Ninit = F.alpha / v0 * F.V\
	  * 1.0 / s0 / sqrt(2.0 * pi)\
	  * exp(- (v - v0) ** 2 / 2 / s0 ** 2)

      pbe_solutions[0] = SteadyStateSolution(
	  Ninit, t, dv,
	  Q=F.Q,
	  gamma=F.gamma,
	  beta=F.beta,
	  pdf='number'
      )

      N = pbe_solutions[0].solution
      m1 = sum(N[:] * v[:])
      m1Init = sum(Ninit[:] * v[:])
      norm = sum(N)
      meanV = m1 / norm
      dMean = (6.0 * meanV / pi) ** (1.0 / 3.0)
      print 'calculated, experimental d: ', dMean, F.expectedD
      print 'mass conservations: ', m1 / m1Init
      Er += sqrt((dMean - F.expectedD) ** 2) / F.expectedD
      #fig = plt.figure()
      #ax = fig.gca()
      #ax.plot(
	  #v, N, "+",
	  #label="result")
      #ax.plot(
	  #v, Ninit, "+",
	  #label="Ninit")
    return Er

#C = fmin(case_error, C0, full_output=True, maxiter=2)

# this was the first guess
#C0 = np.array([10.31e-02, 0.15, 1.0e-12, 5.0e+13])
#----------------simmonsAzzopardi
# after first ten iterations:
C0 = np.array([0.1013e-02, 0.1507, 1.05e-12, 4.78e+13])
# after another 100 iterations:
C0 = np.array([  1.52138957e-03,   6.78265926e-02,   1.05895753e-12,
	      5.12743336e+13])
# -----------------------------------------------------------------

#----------------coulaloglou 0 0
# after first ten iterations (error=0.1):
C0 = np.array([  1.97857877e-03,   6.43129636e-02,   5.27400842e-14,
           1.30520127e+14])
# after another 20 iterations (error = 0.048):
C0 = np.array([  2.51610210e-03,   2.84504679e-02,   5.45292101e-14,
           1.31327842e+14])
# another 20, error 0.0024
C0 = np.array([  2.97439339e-03,   2.54818027e-02,   6.01552014e-14,
           1.02488161e+14])
# another 10, error = 0.0008
C0 = np.array([  2.98817729e-03,   2.55546776e-02,   6.19345476e-14,
           1.02703121e+14])

#----------------coulaloglou 0 1
# after first ten iterations (error=0.029):
C0 = np.array([  3.32158412e-03,   2.60820236e-02,   4.53370793e-14,
           1.15320673e+14])
# after 20 iterations (error=0.0011):
C0 = np.array([  2.81813200e-03,   2.80142045e-02,   4.06690335e-14,
           1.32518719e+14])
# another 20, error = 0.00026
C0 = np.array([  2.84090572e-03,   2.81571210e-02,   4.09260122e-14,
           1.33677217e+14])

#----------------coulaloglou 0 1
# after first 20 iterations (error=0.0039):
C0 = np.array([  2.94435465e-03,   3.13751284e-02,   2.59375997e-14,
           1.70021676e+14])
# after another 20 iterations (error=0.000245):
C0 = np.array([  2.92874399e-03,   3.15027521e-02,   2.67585407e-14,
           1.71109913e+14])

C = fmin(case_error, C0, full_output=True, maxiter=50)
#C = case_error(C0)
#plt.show()
print C
