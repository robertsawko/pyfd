from numpy import arange, zeros
from scipy.integrate import odeint

"""
Method of classes
"""


class MOCSolution:
    """
    Based on Brooks and Hidy uniform discretisation

    """
    def RHS(
        self, N, t
    ):
        dim = N.shape[0]
        # Death breakup term
        dNdt = - N * self.gamma(self.xi)
        # Birth breakup term
        for i in arange(dim):
            for j in arange(i + 1, dim):
                dNdt[i] += \
                    self.beta(self.xi[i], self.xi[j]) \
                    * self.gamma(self.xi[j]) \
                    * N[j] * self.delta_xi

        if self.Q is not None:
            # Birth coalescence term
            for i in arange(dim):
                for j in arange(0, i):
                    dNdt[i] += 0.5 * N[i - j] * N[j] \
                        * self.Q(self.xi[j], self.xi[i - j])
            # Death coalescence term
            for i in arange(dim):
                for j in arange(dim):
                    dNdt[i] -= N[i] * N[j] * self.Q(self.xi[j], self.xi[j])
        return dNdt

    def __init__(
            self, number_of_classes, t, xi0,
            beta=lambda x, y: 2.0 / y, gamma=lambda x: x**2, Q=None
    ):
        # Initial conditions
        N0 = zeros(number_of_classes)
        N0[-1] = 1
        # Kernels setup
        self.beta = beta  # Daughter particle distribution
        self.gamma = gamma  # Breakup frequency
        self.Q = Q  #
        # Uniform grid
        self.xi = xi0 + xi0 * arange(number_of_classes)
        self.delta_xi = xi0
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
