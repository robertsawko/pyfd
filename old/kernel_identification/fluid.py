import numpy as np
from numpy import arange, sum, exp, linspace, sqrt, pi, zeros
import sys
import configparser as ConfigParser


class fluid:
    """
        class to store properties needed for breakup and coalescence models

    """

    def __init__(self, name, caseNr=0, caseNr2=0):
        self.name = name
        self.index = [caseNr, caseNr2]
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
                * self.Re / self.R
            self.Ca = self.muc * self.U[caseNr] / self.sigma
        elif name == "simmonsAzzopardi":
            self.caseNs = np.array([1])
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
            self.epsilon = 0.07978
            self.theta = None
            self.Re = 85000.0
            self.timeRange = arange(0.0, 300.0, 1e-01)
            self.expectedD = 0.32e-03
            #self.d0 = np.array([0.9 * self.expectedD, 1.1 * self.expectedD])
            self.d0 = np.array([0.8 * self.expectedD, 1.2 * self.expectedD])
            self.v0 = pi / 6.0 * self.d0 ** 3
            self.s0 = self.v0 / 4.0
            self.vMax = self.v0 * 6.0
            self.numberOfClasses = 80
            self.R = 2.0 * self.rhoc / (2.0 * self.rhod + self.rhoc)
            self.St = 2.0 / 9.0 * (self.expectedD / self.D) ** 2 \
                * self.Re / self.R
            self.Ca = self.muc * self.U / self.sigma
        elif name == "angeli":
            self.caseNs = np.array([5])
            self.rhoc = 801.0
            self.muc = 1.6e-03
            self.rhod = 1000.0
            self.mud = 1.0e-03
            self.sigma = 1.7e-02
            self.Us = [1.1, 1.25, 1.3, 1.5, 1.7]
            self.U = self.Us[caseNr]
            # volume fraction
            self.alphas = [9.1, 5.0, 7.7, 3.4, 5.0]
            self.alpha = self.alphas[caseNr] * 1e-02
            # pipe diameter, length and volume
            self.D = 24.0e-03
            self.L = 4.5
            self.V = pi * (self.D / 2.0) ** 2 * self.L
            self.epsilons = [0.0075, 0.01, 0.0116, 0.01686, 0.0234] 
            self.epsilon = self.epsilons[caseNr] 
            self.theta = None
            self.Res = [34693.0, 39424.0, 41000.0, 47300.0, 53600.0]
            self.Re = self.Res[caseNr]
            self.timeRange = arange(0.0, 300.0, 1e-01)
            self.expectedDs = [1.29, 0.97, 1.25, 0.73, 0.69]
            self.expectedD = self.expectedDs[caseNr] * 1e-03
            self.d0 = np.array([0.8 * self.expectedD, 1.2 * self.expectedD])
            self.v0 = pi / 6.0 * self.d0 ** 3
            self.s0 = self.v0 / 4.0
            self.vMax = self.v0 * 4.0
            self.numberOfClasses = 50
            self.R = 2.0 * self.rhoc / (2.0 * self.rhod + self.rhoc)
            self.St = 2.0 / 9.0 * (self.expectedD / self.D) ** 2 \
                * self.Re / self.R
            self.Ca = self.muc * self.U / self.sigma
        elif name == "karabelas":
            self.caseNs = np.array([4, 3])
            self.runnumber = []
            self.runnumber.append(['3S', '5S', '6S', '7S'])
            self.runnumber.append(['11S', '17S', '10S'])
            self.rhocs = [798.0, 890.0]
            self.rhoc = self.rhocs[caseNr]
            self.mucs = [1.8e-03, 16.0e-03]
            self.muc = self.mucs[caseNr]
            self.rhod = 1000.0
            self.mud = 1.0e-03
            self.sigmas = [33.1e-03, 34.0e-03]
            self.sigma = self.sigmas[caseNr]
            self.Us = []
            self.Us.append([2.98, 2.57, 2.22, 1.84])
            self.Us.append([3.00, 2.60, 2.24])
            self.U = self.Us[caseNr][caseNr2]
            # volume fraction
            self.alpha = 0.2
            # pipe diameter, length and volume
            self.D = 5.04e-02
            self.L = 4.5
            self.V = pi * (self.D / 2.0) ** 2 * self.L
            self.epsilons = []
            self.epsilons.append([0.145, 0.09855, 0.0671, 0.041])
            self.epsilons.append([0.322, 0.221, 0.150])
            self.epsilon = self.epsilons[caseNr][caseNr2]
            self.theta = None
            self.Res = []
            self.Res.append([66585.0, 57424.0, 49603.0, 41112.0])
            self.Res.append([8410.5, 7289.0, 6280.0])
            self.Re = self.Res[caseNr][caseNr2]
            self.timeRange = arange(0.0, 300.0, 1e-01)
            self.expectedDs = []
            self.expectedDs.append([407.0, 455.0, 794.0, 1066.0])
            self.expectedDs.append([347.0, 401.0, 512.0])
            self.expectedD = self.expectedDs[caseNr][caseNr2] * 1e-06
            self.d0 = np.array([0.8 * self.expectedD, 1.2 * self.expectedD])
            self.v0 = pi / 6.0 * self.d0 ** 3
            self.s0 = self.v0 / 4.0
            self.vMax = self.v0 * 4.0
            self.numberOfClasses = 80
            self.R = 2.0 * self.rhoc / (2.0 * self.rhod + self.rhoc)
            self.St = 2.0 / 9.0 * (self.expectedD / self.D) ** 2 \
                * self.Re / self.R
            self.Ca = self.muc * self.U / self.sigma
        elif name == "coulaloglou":
            self.caseNs = np.array([5, 5, 4])
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
                * self.Re / self.R
            self.Ca = self.mud * self.Nstar * self.Dstar / self.sigma \
                * sqrt(self.rhoc / self.rhod)
        else:
            sys.exit("Valid cases are: 'galinat', 'simmonsAzzopardi', 'coulaloglou', 'angeli', 'karabelas'")

        # default values from Coulaloglou and Tavlarides
        self.C = np.array([0.04, 0.08, 2.17e-16, 2.28e13])

    def gamma(self, xi):
        C = self.C[0] * xi ** (-2.0 / 9.0) * self.epsilon ** (1.0 / 3.0)\
            / (1.0 + self.alpha)
            #* self.Re / (1.0 + self.alpha)

        exp_argument = - self.C[1] * self.sigma * (1.0 + self.alpha) ** 2 \
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

        exp_argument = - self.C[3] * self.muc * self.rhoc * self.epsilon\
            / (1.0 + self.alpha) ** 3 * dRatio / self.sigma ** 2

        C = self.C[2] * (xi1 ** (1.0 / 3.0) + xi2 ** (1.0 / 3.0)) ** 2\
            * (xi1 ** (2.0 / 9.0) + xi2 ** (2.0 / 9.0)) ** 0.5\
            * self.epsilon ** (1.0 / 3.0) / (1.0 + self.alpha) / self.V
        return exp(exp_argument) * C

    def We(self):
        return 2.0 * self.epsilon ** (2.0 / 3.0) * self.rhoc\
            * self.expectedD ** (5.0 / 3.0) / self.sigma

    def initialC(self):
        data = ConfigParser()
        data.read('validationData/initial_constants.ini')
        name = self.name + '-'\
            + repr(self.index[0]) + '-' + repr(self.index[1])
        C = []
        for i in range(4):
            C.append((float)(data.get(name, 'C' + repr(i + 1))))
        return C
