from numpy import arange, sum, exp, sqrt, pi
from steadyStateSolver import SteadyStateSolution
from fluid import fluid
import numpy as np
from scipy.optimize import fmin
from sys import argv

Er = 0
error_file = "validationData/convergence_example.txt"


def write_error(C):
    with open(error_file, "a") as myfile:
            myfile.write(repr(Er))
            for c in C:
                myfile.write(" " + repr(c))
            myfile.write('\n')


def case_error(C):
    F = fluid(argv[1], (int)(argv[2]), (int)(argv[3]))

    F.C = C
    print '--------------'
    print 'constants: ', F.C[0], F.C[1], F.C[2], F.C[3]

    global Er
    Er = 0.0
    for i in range(2):
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

        pbe_solutions = SteadyStateSolution(
            Ninit, t, dv,
            Q=F.Q,
            gamma=F.gamma,
            beta=F.beta,
            pdf='number'
        )
        pbe_solutions.solve()

        N = pbe_solutions.solution
        m1 = sum(N[:] * v[:])
        m1Init = sum(Ninit[:] * v[:])
        norm = sum(N)
        meanV = m1 / norm
        dMean = (6.0 * meanV / pi) ** (1.0 / 3.0)
        print 'calculated, experimental d: ', dMean, F.expectedD
        print 'mass conservations: ', m1 / m1Init
        Er += sqrt((dMean - F.expectedD) ** 2) / F.expectedD
    return Er

# -----------------------------------------------------------------


# TODO: move C0 to fluid class
F = fluid(argv[1], (int)(argv[2]), (int)(argv[3]))
C0 = np.array(F.initialC())
#C0 = 0.75 * C0

# TODO: construct list of cases
open(error_file, 'w').close()
C = fmin(case_error, C0, full_output=True, maxiter=50, callback=write_error)
print C
# TODO write a for loop over the list of cases with parallel execution
