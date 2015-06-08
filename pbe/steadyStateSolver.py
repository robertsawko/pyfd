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


    def solve(self, N):
      rhs = self.RHS(N, 0.0)
      er0 = sum(rhs * rhs)
      N0 = N

      i = 0
      iMax = 10
      iMin = 10
      er = er0
      #for i in arange(rhs.shape[0]):
	#print rhs[i], N0[i]
      m1Init = sum(self.xi * N0)
      #print m1Init

      t = np.arange(0, 5.0, 1.0)
      while er > er0 * 1e-03 and i < iMax or i < iMin:
	i += 1
	Ni = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
	N0 = Ni[-1]
	m1 = sum(self.xi * N0)
	x = m1 / m1Init
	if x > 1.4 or x < 0.6:
	  return N0
	rhs = self.rhs
	er = sum(rhs * rhs)
	#print er
      print "Converged in ", i , " iterations"
      return N0


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
        # Uniform grid
        self.xi = xi0 + xi0 * arange(self.number_of_classes)
        self.delta_xi = xi0
	self.solution = self.solve(N0)
	self.dt = 0.1
        # Solve procedure
        #self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
	#tol = sum(N0) / N0.shape[0] * 0.1
	#print tol
	#self.solution = broyden1(self.RHS, N0, verbose=True, f_tol=tol)
	#print.solution = fsolve(self.RHS, N0, maxfev=2, full_output=True)
