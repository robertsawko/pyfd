import numpy as np
from numpy import arange, sum, exp, linspace, sqrt, pi, zeros
import sys


class fluid:
    """
        class to store properties needed for breakup and coalescence models

    """

    def __init__(self, name, caseNr=0, caseNr2=0):
        if name == "galinat":
            self.rhoc = 996.0
            self.rhod = 683.7
            self.mud = 4.5e-04
            self.muc = 8.2e-04
            self.sigma = 4.7e-02
            self.dpMax = np.array([173.0, 363.0, 566.0, 706.0, 871.0, 1120.0])
            self.U = np.array([0.118, 0.177, 0.236, 0.275, 0.314, 0.354])
            self.Res = np.array([4250, 6400, 8600, 10000, 11500, 12900])
            # orifice ratio
            self.beta_or = 0.5
            # volume fraction
            self.alpha = 0.02
            # pipe diameter, and length
            self.D = 0.03
            self.L = 1.0
            self.Re = self.Res[caseNr]
            # residence time; this is a though one...
            # we'll take it equal to the time needed for the
            # mean flow to travel length equal to three time
            # the orifice thickness
            # orifice thickness is 5mm
            self.thetas = 2.0 * 0.005 / self.U[:]
            self.epsilons = 1.0 / self.rhoc * self.dpMax[:] * self.U[:]\
                / 2.0 / self.D * (1.0 / self.beta_or ** 2 - 1.0)
            self.epsilon = self.epsilons[caseNr]
            self.theta = self.thetas[caseNr]
            self.V = pi * (self.D / 2.0) ** 2 * self.L
            self.timeRange = arange(0.0, 10.0 * self.theta, 1e-02)
            self.d0s = np.genfromtxt('validationData/orifice/d32_upstream.txt',
                                     delimiter=',')
            self.d0 = self.d0s[caseNr][1] * 1e-03
            exp_dRe =\
                np.genfromtxt('validationData/orifice/d32_downstream.txt',
                              delimiter=',')
            self.expectedD = exp_dRe.T[1][caseNr] * 1e-03
            self.v0 = pi / 6.0 * self.d0 ** 3
            self.s0 = self.v0 / 8.0
            self.vMax = self.v0 * 1.3
            self.numberOfClasses = 120
            self.R = 2.0 * self.rhoc / (2.0 * self.rhod + self.rhoc)
            self.St = 2.0 / 9.0 * (self.expectedD / self.D) ** 2 \
                * self.Re * self.R
            self.Ca = self.muc * self.U[caseNr] / self.sigma
        elif name == "simmonsAzzopardi":
            self.rhoc = 797.0
            self.muc = 1.8e-03
            self.rhod = 1166.0
            self.mud = 1.6e-03
            self.sigma = 1.0e-02
            self.U = 2.71
            # volume fraction
            self.alpha = 0.117
            # pipe diameter, length and volume
            self.D = 0.063
            self.L = 4.5
            self.V = pi * (self.D / 2.0) ** 2 * self.L
            self.epsilon = 0.082
            self.theta = None
            self.Re = 78200.0
            self.timeRange = arange(0.0, 300.0, 1e-01)
            self.expectedD = 0.32e-03
            #self.d0 = np.array([0.9 * self.expectedD, 1.1 * self.expectedD])
            self.d0 = np.array([0.8 * self.expectedD, 1.2 * self.expectedD])
            self.v0 = pi / 6.0 * self.d0 ** 3
            self.s0 = self.v0 / 4.0
            self.vMax = self.v0 * 4.0
            self.numberOfClasses = 50
            self.R = 2.0 * self.rhoc / (2.0 * self.rhod + self.rhoc)
            self.St = 2.0 / 9.0 * (self.expectedD / self.D) ** 2 \
                * self.Re * self.R
            self.Ca = self.muc * self.U / self.sigma
        elif name == "coulaloglou":
            self.rhoc = 1000.0
            self.muc = 1.0e-03
            self.rhod = 972.0
            self.mud = 1.3e-03
            self.sigma = 42.82e-03
            self.alphas = np.array([0.05, 0.1, 0.15])
            self.alpha = self.alphas[caseNr]
            exp_dRe =\
                np.genfromtxt('validationData/coulaloglou/d32_N_alpha'
                              + repr((caseNr + 1) * 5) + '.txt')
            self.expectedD = exp_dRe.T[1][caseNr2] * 1e-03
            self.Nstar = exp_dRe.T[0][caseNr2] / 60.0
            # impeller speed: (in 1 / min)
            #self.Nstars = np.array([190.0, 220.0, 250.0, 280.0, 310.0])
            # convert to 1 / second:
            #self.Nstar = self.Nstars[caseNr2] / 60.0
            # impeller diameter: (10cm)
            self.Dstar = 0.1
            # tank volume (12l)
            self.V = 12.0e-03
            self.epsilon = 0.407 * self.Nstar ** 3 * self.Dstar ** 2
            # residence time is 10 minutes
            self.theta = 10.0 * 60.0
            self.Re = self.Nstar * self.Dstar ** 2 / self.muc * self.rhoc
            self.timeRange = arange(0.0, 1.0, 1e-02)
            self.d0 = np.array([0.9 * self.expectedD, 1.1 * self.expectedD])
            #self.expectedD = 0.255e-03
            self.v0 = pi / 6.0 * self.d0 ** 3
            self.s0 = self.v0 / 3.0
            self.vMax = self.v0 * 12.0
            self.numberOfClasses = 70
            self.R = 2.0 * self.rhoc / (2.0 * self.rhod + self.rhoc)
            self.St = 2.0 / 9.0 * (self.expectedD / self.Dstar) ** 2 \
                * self.Re * self.R
            self.Ca = self.mud * self.Nstar * self.Dstar / self.sigma \
                * sqrt(self.rhoc / self.rhoc)
        else:
            sys.exit("Valid cases are: 'galinat', 'simmonsAzzopardi', 'coulaloglou'")

        # default values from Coulaloglou and Tavlarides
        self.C1 = 0.04
        self.C2 = 0.08
        self.C3 = 2.17e-16
        self.C4 = 2.28e13

    def gamma(self, xi):
        C = self.C1 * xi ** (-2.0 / 9.0) * self.epsilon ** (1.0 / 3.0)\
            / (1.0 + self.alpha)
            #* self.Re / (1.0 + self.alpha)

        exp_argument = - self.C2 * self.sigma * (1.0 + self.alpha) ** 2 \
            / (self.rhod * xi ** (5.0 / 9.0) * self.epsilon ** (2.0 / 3.0))
            #/ self.Re
        return C * exp(exp_argument)

    # droplet daughter distribution:
    def beta(self, xi2, xi1):
        return 2.0 * 2.4 / xi1 * exp(- 4.5 * (2.0 * xi2 - xi1) ** 2 / xi1 ** 2)

    # coalescence rate:
    def Q(self, xi1, xi2):
        dRatio = xi1 ** (1.0 / 3.0) * xi2 ** (1.0 / 3.0)\
            / (xi1 ** (1.0 / 3.0) + xi2 ** (1.0 / 3.0))
        dRatio = dRatio ** 4

        exp_argument = - self.C4 * self.muc * self.rhoc * self.epsilon\
            / (1.0 + self.alpha) ** 3 * dRatio / self.sigma ** 2

        C = self.C3 * (xi1 ** (1.0 / 3.0) + xi2 ** (1.0 / 3.0)) ** 2\
            * (xi1 ** (2.0 / 9.0) + xi2 ** (2.0 / 9.0)) ** 0.5\
            * self.epsilon ** (1.0 / 3.0) / (1.0 + self.alpha) / self.V
        return exp(exp_argument) * C
