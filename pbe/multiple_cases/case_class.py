from moc import MOCSolution
from numpy import arange, sqrt, exp, pi, zeros


def beta(v1, v2):
    return 2.4 / v2 * exp(-4.5 * (2 * v1 - v2)**2 / (v2**2))


class CaseSolution(MOCSolution):
    def __init__(
            self,
            dispProperties,
            contProperties,
            domainProperties,
            modelParameters=None):
        self.dispProperties = dispProperties
        self.contProperties = contProperties
        self.phi = dispProperties['phi']

        self.muc = contProperties['mu']
        self.rhoc = contProperties['rho']
        self.epsilon = contProperties['epsilon']
        self.rhod = dispProperties['rho']
        self.sigma = dispProperties['sigma']

        time = arange(0.0, 3600, 1)

        vmax = dispProperties['vMax']

        # Feed distribution
        self.v0 = dispProperties['v0']
        self.sigma0 = dispProperties['sigma0']

        # Feed
        theta = domainProperties['theta']
        M = domainProperties['M']
        self.Vt = domainProperties['V']
        if theta is None:
            self.n0 = None
        else:
            self.n0 = self.Vt / theta

        Ninit = zeros(M)

        if modelParameters is None:
            self.C1 = 0.4
            self.C2 = 0.08
            self.C3 = 2.8e-6
            self.C4 = 1.83e9
        else:
            self.C1 = modelParameters['C1']
            self.C2 = modelParameters['C2']
            self.C3 = modelParameters['C3']
            self.C4 = modelParameters['C4']

        MOCSolution.__init__(
            self, Ninit, time, vmax / M,
            beta=beta, gamma=self.g, Q=self.Q,
            theta=theta, n0=self.n0, A0=self.A0)

    def A0(self, v):
        return \
            self.phi / (self.v0 * self.sigma0 * sqrt(2 * pi)) * \
            exp(-(v - self.v0)**2 / (2 * self.sigma0**2))

    def g(self, v):
        C1 = self.C1
        C2 = self.C2
        return \
            C1 * v**(-2. / 9) * self.epsilon**(1. / 3) / (1 + self.phi) * \
            exp(- C2 * (1 + self.phi)**2 * self.sigma /
                (self.rhod * v**(5. / 9) * self.epsilon**(2. / 3)))

    def Q(self, v1, v2):
        C3 = self.C3
        C4 = self.C4
        d_ratio = (v1**(1. / 3) * v2**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))

        return C3 * \
            (v1**(1. / 3) + v2**(1. / 3))**2. * \
            (v1**(2. / 9) + v2**(2. / 9))**0.5 * \
            self.epsilon**(1. / 3) / \
            ((1 + self.phi) * self.Vt) * \
            exp(
                -C4 * self.muc * self.rhoc * self.epsilon /
                self.sigma**2 /
                (1 + self.phi)**3 *
                d_ratio**4)

    @property
    def pbe_phi(self):
        return self.total_volume / self.Vt
