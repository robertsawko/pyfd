from numpy import arange, zeros
from scipy.integrate import odeint
import sys
from scipy.optimize import newton_krylov, fsolve, broyden1
import numpy as np

"""
Method of classes
"""


class SteadyStateSolution:
    """
    Based on Brooks and Hidy uniform discretisation

    """
    def RHS(
        self, N, t
    ):
       #pdb.set_trace()
        dNdt = zeros(self.number_of_classes)

        # TODO: Refactor that to MOC and create MOC inheritance
        if self.gamma is not None and self.beta is not None:
            for i in arange(self.number_of_classes):
                # Death breakup term
                if i != 0:
                    dNdt[i] -= N[i] * self.gamma(self.xi[i])
                # Birth breakup term
                if i != (self.number_of_classes - 1):
                    for j in arange(i + 1, self.number_of_classes):
                        dNdt[i] += \
                            self.beta(self.xi[i], self.xi[j]) \
                            * self.gamma(self.xi[j]) \
                            * N[j] * self.delta_xi

        if self.Q is not None:
            for i in arange(self.number_of_classes):
                # Birth coalescence term
                if i != 0:
                    for j in arange(i):
                        dNdt[i] += 0.5 * N[i - j - 1] * N[j] \
                            * self.Q(self.xi[j], self.xi[i - j - 1])
                # Death coalescence term
                if i != (self.number_of_classes - 1):
                    for j in arange(self.number_of_classes):
                        dNdt[i] -= N[i] * N[j] * self.Q(self.xi[i], self.xi[j])

        if self.theta is not None:
            dNdt -= (N - self.N0) / self.theta
        self.rhs = dNdt
        return dNdt

    def setIter(self, iMin, iMax):
        self.iMin = iMin
        self.iMax = iMax

    def setMassTol(self, mMin, mMax):
        self.mMin = mMin
        self.mMax = mMax

    def setTime(self, t):
        self.time = t

    def setError(self, er):
        self.er = er

    def solve(self):
        N0 = self.N0
        rhs = self.RHS(N0, 0.0)
        er0 = sum(rhs * rhs)

        i = 0
        iMax = self.iMax
        iMin = self.iMin
        er = er0
        m1Init = sum(self.xi * N0)

        t = self.time
        while er > er0 * self.er and i < iMax or i < iMin:
            i += 1
            Ni = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
            N0 = Ni[-1]
            m1 = sum(self.xi * N0)
            x = m1 / m1Init
            if x > self.mMax or x < self.mMin:
                break
            rhs = self.rhs
            er = sum(rhs * rhs)
        print("Converged in ", i, " iterations")
        #return N0
        self.solution = N0

    def __init__(self, N0, t, xi0, beta=None, gamma=None, Q=None, theta=None,
                 pdf="number"):
        self.number_of_classes = N0.shape[0]
        # Kernels setup
        self.beta = beta  # Daughter particle distribution
        self.gamma = gamma  # Breakup frequency
        self.Q = Q  #
        # inflow and outflow replaced with relaxation to equilibrium
        # process with relaxation time equal to residence time theta
        self.theta = theta
        self.N0 = N0
        # choose pdf formulation: number or number density
        self.pdf = pdf
        # set minimum and maximum number of iteration:
        self.setIter(10, 100)
        self.setMassTol(0.6, 1.4)
        self.setTime(np.arange(0, 0.5, 0.1))
        self.setError(1e-03)
        # Uniform grid
        self.xi = xi0 + xi0 * arange(self.number_of_classes)
        self.delta_xi = xi0
        self.N0 = N0
        #self.solution = self.solve(N0)
