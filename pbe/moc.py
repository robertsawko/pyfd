from numpy import arange, zeros, pi, zeros_like, dot, array
from numpy import sum as nsum
from scipy.integrate import odeint, quad

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
        dNdt = zeros_like(N)

        if self.gamma is not None and self.betadxi is not None:
            # Death breakup term
            dNdt[1:] -= N[1:] * self.gamma[1:]
            dNdt[:-1] += self.nu * dot(
                self.betadxi[:-1, 1:], N[1:] * self.gamma[1:])

        if self.Q is not None:
            Cd = zeros_like(dNdt)
            for i in arange(self.number_of_classes / 2):
                ind = slice(i, self.number_of_classes - i - 1)
                Cb = self.Q[i, ind] * N[i] * N[ind]
                Cd[i] += nsum(Cb)
                Cd[(i + 1):(i + len(Cb))] += Cb[1:]
                Cb[0] = 0.5 * Cb[0]
                dNdt[(2 * i + 1):] += Cb

            dNdt -= Cd
        if self.theta is not None:
            dNdt += (self.n0 * self.A0 - N / self.theta)

        # print('Time = {0:g}'.format(t))
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
            N0 = zeros_like(self.xi)
        else:
            N0 = array([
                quad(N0, self.xi[i] - dxi / 2., self.xi[i] + dxi / 2.)[0]
                for i in range(number_of_classes)])

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
                    self.betadxi[:, i] / nsum(self.betadxi[:, i])

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

        if A0 is None:
            self.A0 = None
        else:
            self.A0 = array([
                quad(A0, self.xi[i] - dxi / 2., self.xi[i] + dxi / 2.)[0]
                for i in range(number_of_classes)])
        # Solve procedure
        self.N = odeint(lambda NN, t: self.RHS(NN, t), N0, t)
