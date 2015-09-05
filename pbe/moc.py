from numpy import arange, zeros, pi
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
            dNdt[1:] -= N[1:] * self.gamma[1:]
            for i in arange(self.number_of_classes):
                # Birth breakup term
                if i != (self.number_of_classes - 1):
                    for j in arange(i + 1, self.number_of_classes):
                        dNdt[i] += \
                            self.beta[i, j] * \
                            self.gamma[j] * \
                            N[j] * self.delta_xi

        if self.Q is not None:
            for i in arange(self.number_of_classes):
                # Birth coalescence term
                if i != 0:
                    for j in arange(i):
                        dNdt[i] += 0.5 * N[i - j - 1] * N[j] *\
                            self.Q[j, i - j - 1]
                # Death coalescence term
                if i != (self.number_of_classes - 1):
                    for j in arange(self.number_of_classes):
                        dNdt[i] -= N[i] * N[j] * self.Q[i, j]

        if self.theta is not None:
            dNdt += (self.n0 * self.A0 - N / self.theta)

        return dNdt

    @property
    def number_density(self):
        return self.N / self.xi0

    @property
    def d32(self):
        return \
            (6 / pi * sum(self.N[-1] * self.xi) / sum(self.N[-1]))**(1. / 3)

    @property
    def total_volume(self):
        return sum(self.N[-1] * self.xi)

    def __init__(
            self, N0, t, xi0,
            beta=None, gamma=None, Q=None,
            theta=None, n0=None, A0=None):
        self.number_of_classes = N0.shape[0]
        self.N0 = N0
        self.xi0 = xi0
        self.n0 = n0
        self.A0 = A0
        self.theta = theta
        # Uniform grid
        self.xi = xi0 + xi0 * arange(self.number_of_classes)
        # Kernels setup
        if gamma is not None:
            self.gamma = gamma(self.xi)
            self.beta = zeros((self.number_of_classes, self.number_of_classes))
            for i in range(len(self.xi)):
                for j in range(len(self.xi)):
                    self.beta[i, j] = beta(self.xi[i], self.xi[j])

        else:
            self.gamma = None
            self.beta = None

        if Q is not None:
            self.Q = zeros((self.number_of_classes, self.number_of_classes))
            for i in range(len(self.xi)):
                for j in range(len(self.xi)):
                    self.Q[i, j] = Q(self.xi[i], self.xi[j])  #
        else:
            self.Q = None
        self.delta_xi = xi0
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
