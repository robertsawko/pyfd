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
            Re = F.Re
            Ca = F.Ca
            St = F.St
            We = F.We()
            # -----------------------------
            A = 0.11371549951931224
            B = 0.0020814470588332158
            aWe = 1.1466914631499581
            aCa = -0.14907380395843822
            C1 = A * We ** aWe / Ca ** aCa + B
            #C1 = 0.0039507407407407408
            # -----------------------------
            A = 31.563299222959419
            B = 0.019272087314426867
            aWe = 1.6170490368507624
            aRe = -0.2238711406747752
            C2 = A * We ** aWe * Re ** aRe + B
            #C2 = 0.035507407405302893
            # -----------------------------
            A = 1.6428900960606153e-13
            B = 6.1228384797660674e-15
            aCa = 1.3927752887548224
            aSt = 0.71738786775149754
            c1 = 3.5541424951605571
            c2 = 0.63108186158198221
            C3 = A * St ** aSt * (Ca ** aCa * c1 + c2) + B
            # -----------------------------
            A = 513765273.33555937
            B = 0.0
            aWe = 0.18740863508646571
            aSt = -0.6257571490782956
            C4 = A * Re * St ** aSt / We ** aWe

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
