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
            dNdt += (self.u0 * self.N0 - N / self.theta)

        return dNdt

    def number_density(self):
        return self.N / self.xi0

    def __init__(
            self, N0, t, xi0,
            beta=None, gamma=None, Q=None,
            theta=None, u0=None):
        self.number_of_classes = N0.shape[0]
        self.N0 = N0
        # Kernels setup
        self.beta = beta  # Daughter particle distribution
        self.gamma = gamma  # Breakup frequency
        self.Q = Q  #
        self.xi0 = xi0
        self.u0 = u0
        self.theta = theta
        # Uniform grid
        self.xi = xi0 + xi0 * arange(self.number_of_classes)
        self.delta_xi = xi0
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
