import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
from itertools import cycle
import fluid
from steadyStateSolver import SteadyStateSolution

# use 0:
# dinit = 0.8 * d_exp
# use 1:
# dinit = 1.2 * d_exp
n = 0

dNum = []
dExp = []

caseName = ['coulaloglou', 'angeli', 'karabelas', 'simmonsAzzopardi']
caseNs = []

for name in caseName:
  F = fluid.fluid(name, 0, 0)
  caseNs.append(F.caseNs)

caseId = 0

for case in caseNs:
  I = case.shape[0]
  for i in np.arange(I):
    for j in np.arange(case[i]):
      pbe_solutions = dict()
      F = fluid.fluid(caseName[caseId], i, j)
      C1 = -3.15815e-06 * F.Re ** 0.598 / F.St ** 0.1502 / F.Ca ** 0.06395 + 0.0063983
      C2 = -8.46752e-08 * F.Re ** 1.11306 + 0.0473
      C3 = 1.6268e-12 * 0.062 * np.log(F.Re) * F.St ** 0.65776 * F.Ca ** 0.54667 + 1.7367e-14
      C4 = 284255 * F.Re ** 1.6296 / F.St ** 0.5044 / F.Ca ** 0.2987
      F.C = np.array([C1, C2, C3, C4])
      print F.C
      t = F.timeRange
      v0 = F.v0[n]
      s0 = F.s0[n]
      vmax = F.vMax[n]
      g = F.numberOfClasses
      print "initial d and v: ", i, F.d0[n], v0

      dv = vmax / g
      v = dv + dv * np.arange(g)
      Ninit = F.alpha / v0 * F.V\
	  * 1.0 / s0 / np.sqrt(2.0 * np.pi)\
	  * np.exp(- (v - v0) ** 2 / 2 / s0 ** 2)

      pbe_solutions[0] = SteadyStateSolution(
	  Ninit, t, dv,
	  Q=F.Q,
	  gamma=F.gamma,
	  beta=F.beta,
	  pdf='number'
      )
      pbe_solutions[0].solve()

      N = pbe_solutions[0].solution
      m1 = sum(N[:] * v[:])
      m1Init = sum(Ninit[:] * v[:])
      norm = sum(N)
      meanV = m1 / norm
      dMean = (6.0 * meanV / np.pi) ** (1.0 / 3.0)
      print 'calculated, experimental d: ', dMean, F.expectedD
      print 'mass conservations: ', m1 / m1Init
      print '------------------------'
      dNum.append(dMean)
      dExp.append(F.expectedD)

  #------------------------------------------------
  caseId += 1

np.savetxt('validationData/comparison/dExp.txt', dExp)
np.savetxt('validationData/comparison/dNum.txt', dNum)
