from .moc import MOCSolution
from numpy import arange, sqrt, exp, pi


class DomainProperties:
    def __init__(self, theta, V, M):
        self.theta = theta
        self.V = V
        self.M = M


class DispersionProperties:
    def __init__(self, phi, rho, sigma, v_max, v0, sigma0):
        self.phi = phi
        self.rho = rho
        self.sigma = sigma
        self.v_max = v_max
        self.v0 = v0
        self.sigma0 = sigma0


def beta(v1, v2):
    return 2.4 / v2 * exp(-4.5 * (2 * v1 - v2)**2 / (v2**2))


'''
input:
    dispersion:
        description: is a dictionary with properties of the dispersed phase
        dispersion = dict()
        Required fields:
        phi         volume fraction
        rho         density
        sigma       interfacial tension
        vMax        maximum volume for the discretization
        v0          mean volume for the distributinon function
        sigma0      standard deviation for the distribution function

    contProperties:
        description: is a dictionary with properties of the continuous phase
        Required fields:
        mu          viscosity
        rho         density
        epsilon     turbulent energy dissipation rate

    domainProperties:
        description: is a dictionary with properties of the computational
        domain and discretization parameters
        Required fields:
        theta       mean residence time
        M           number of classes used
        V           domain volume

    model_parameters:
        description: list of four values C1 - C4 representing breakup model
        constants (C1 and C2) and coalescence model constants (C3 and C4)

    time:
        description: discretized time domain
'''


class CaseSolution(MOCSolution):
    def __init__(
            self,
            dispersion,
            contProperties,
            domain,
            model_parameters=None,
            time=arange(0.0, 3600, 0.5)):

        self.contProperties = contProperties
        self.phi = dispersion.phi

        self.muc = contProperties['mu']
        self.rhoc = contProperties['rho']
        self.epsilon = contProperties['epsilon']
        self.rhod = dispersion.rho
        self.sigma = dispersion.sigma

        vmax = dispersion.v_max

        # Feed distribution
        self.v0 = dispersion.v0
        self.sigma0 = dispersion.sigma0

        # Feed
        theta = domain.theta
        M = domain.M  # TODO: Fix variable name, what is M?
        self.Vt = domain.Vt
        if theta is None:
            self.n0 = None
        else:
            self.n0 = self.Vt / theta

        if model_parameters is None:
            self.C = [0.4, 0.08, 2.8, 1.83e13]
        else:
            self.C = model_parameters

        MOCSolution.__init__(
            self, M, time, vmax / M,
            beta=beta, gamma=self.g, Q=self.Qf,
            theta=theta, n0=self.n0, A0=self.A0)

    def A0(self, v):
        return \
            self.phi / (self.v0 * self.sigma0 * sqrt(2 * pi)) * \
            exp(-(v - self.v0)**2 / (2 * self.sigma0**2))

    def g(self, v):
        C1 = self.C[0]
        C2 = self.C[1]
        return \
            C1 * v**(-2. / 9) * self.epsilon**(1. / 3) / (1 + self.phi) * \
            exp(- C2 * (1 + self.phi)**2 * self.sigma /
                (self.rhod * v**(5. / 9) * self.epsilon**(2. / 3)))

    def Qf(self, v1, v2):
        C3 = self.C[2]
        C4 = self.C[3]
        d_ratio = (v1**(1. / 3) * v2**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))

        return C3 * \
            (v1**(2. / 3) + v2**(2. / 3)) * \
            (v1**(2. / 9) + v2**(2. / 9))**0.5 * \
            self.epsilon**(1. / 3) / \
            ((1 + self.phi)) * \
            exp(
                -C4 * self.muc * self.rhoc * self.epsilon /
                self.sigma**2 /
                (1 + self.phi)**3 *
                d_ratio**4)

    @property
    def pbe_phi(self):
        return self.total_volume / self.Vt
