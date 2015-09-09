from numpy import arange, zeros, pi
from numpy import sum as nsum
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
        print("Time = {0:g}".format(t))
        dNdt = zeros(self.number_of_classes)

        if self.gamma is not None and self.betadxi is not None:
            # Death breakup term
            for i in arange(1, self.number_of_classes):
                dNdt[i] -= N[i] * self.gamma[i]
                for j in arange(i):
                    dNdt[j] += \
                        self.nu * self.betadxi[j, i] * \
                        self.gamma[i] * \
                        N[i]

        if self.Q is not None:
            for i in arange(self.number_of_classes):
                # Birth coalescence term
                if i != 0:
                    for j in arange(i):
                        dNdt[i] += 0.5 * N[i - j - 1] * N[j] *\
                            self.Q[j, i - j - 1]
                # Death coalescence term
                if i != (self.number_of_classes - 1):
                    for j in arange(self.number_of_classes - i - 1):
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
    def xi_d(self):
        return \
            (6 / pi * self.xi)**(1. / 3)

    @property
    def total_volume(self):
        return nsum(self.N[-1] * self.xi)

    @property
    def total_numbers(self):
        return nsum(self.N, axis=1)

    def __init__(
            self, number_of_classes, t, dxi, N0=None, xi0=None,
            beta=None, gamma=None, Q=None,
            theta=None, n0=None, A0=None):
        self.number_of_classes = number_of_classes
        if xi0 is None:
            self.xi0 = dxi
        else:
            self.xi0 = xi0
        self.n0 = n0
        self.theta = theta
        # Uniform grid
        self.xi = self.xi0 + dxi * arange(self.number_of_classes)
        if N0 is None:
            N0 = zeros(number_of_classes)
        else:
            N0 = N0(self.xi) * dxi

        self.nu = 2.0  # Binary breakup
        # Kernels setup
        if gamma is not None:
            self.gamma = gamma(self.xi)
            self.betadxi = zeros(
                (self.number_of_classes, self.number_of_classes))
            for i in range(1, len(self.xi)):
                for j in range(i):
                    self.betadxi[j, i] = beta(self.xi[j], self.xi[i])
                self.betadxi[:, i] =\
                    self.betadxi[:, i] / sum(self.betadxi[:, i])

        else:
            self.gamma = None
            self.betadxi = None

        if Q is not None:
            self.Q = zeros((self.number_of_classes, self.number_of_classes))
            for i in range(len(self.xi)):
                for j in range(len(self.xi)):
                    self.Q[i, j] = Q(self.xi[i], self.xi[j])  #
        else:
            self.Q = None

        if A0 is not None:
            self.A0 = A0(self.xi) * dxi
        else:
            self.A0 = None
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
