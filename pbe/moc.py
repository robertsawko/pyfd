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
        self, N, t, xi, deltaXi,
        beta=lambda x, y: 2.0 / y, gamma=lambda x: x**2
    ):
        dim = N.shape[0]
        dNdt = - N * gamma(xi)
        for i in arange(dim):
            for j in arange(i + 1, dim):
                dNdt[i] += beta(xi[i], xi[j]) * gamma(xi[j]) * N[j] * deltaXi
        return dNdt

    def __init__(self, number_of_classes, t, xi0):
        N0 = zeros(number_of_classes)
        N0[-1] = 1
        self.xi = xi0 + xi0 * arange(number_of_classes)
        self.N = odeint(lambda NN, t: self.RHS(NN, t, self.xi, xi0), N0, t)
