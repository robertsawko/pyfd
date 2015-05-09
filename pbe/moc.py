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
        dNdt = zeros(self.number_of_classes)

        if self.gamma is not None and self.beta is not None:
            # Death breakup term
            dNdt -= N * self.gamma(self.xi)
            # Birth breakup term
            for i in arange(self.number_of_classes):
                for j in arange(i + 1, self.number_of_classes):
                    dNdt[i] += \
                        self.beta(self.xi[i], self.xi[j]) \
                        * self.gamma(self.xi[j]) \
                        * N[j] * self.delta_xi

        if self.Q is not None:
            for i in arange(self.number_of_classes):
                # Birth coalescence term
                for j in arange(0, i):
                    dNdt[i] += 0.5 * N[i - j] * N[j] \
                        * self.Q(self.xi[j], self.xi[i - j])
                # Death coalescence term
                for j in arange(self.number_of_classes):
                    dNdt[i] -= N[i] * N[j] * self.Q(self.xi[i], self.xi[j])
        return dNdt

    def __init__(self, N0, t, xi0, beta=None, gamma=None, Q=None):
        self.number_of_classes = N0.shape[0]
        # Kernels setup
        self.beta = beta  # Daughter particle distribution
        self.gamma = gamma  # Breakup frequency
        self.Q = Q  #
        # Uniform grid
        self.xi = xi0 + xi0 * arange(self.number_of_classes)
        self.delta_xi = xi0
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
