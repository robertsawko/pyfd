from numpy import arange, sum, exp, linspace, sqrt, pi, zeros, genfromtxt
from itertools import cycle
from steadyStateSolver import SteadyStateSolution
import matplotlib.pyplot as plt
from fluid import fluid
import os.path
import numpy as np
from scipy.optimize import minimize, fmin, fmin_powell
from sys import argv


def case_error(C):
    pbe_solutions = dict()
    F = fluid(argv[1], (int)(argv[2]), (int)(argv[3]))

    F.C = C
    print '--------------'
    print 'constants: ', F.C[0], F.C[1], F.C[2], F.C[3]

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
        pbe_solutions[0].solve()

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

# -----------------------------------------------------------------

C0 = np.array([6.02352486e-03,   5.44936761e-02,   1.39863073e-13,
             3.35891500e+12])

#C = fmin(case_error, C0, full_output=True, maxiter=50)
#print C
F = fluid('simmonsAzzopardi', 0, 0)
print F.Re, F.St, F.Ca
