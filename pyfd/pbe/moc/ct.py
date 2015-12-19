from .moc import MOCSolution
from numpy import arange, sqrt, exp, pi

"""
Case setup based on:

    Colaloglou and Tavlarides (1977)
    "Description of Interaction Processes in Agitated Liquid-Liquid
    Dispersions", J. Chem Eng Vol 32

"""


def beta(v1, v2):
    return 2.4 / v2 * exp(-4.5 * (2 * v1 - v2)**2 / (v2**2))


class CTSolution(MOCSolution):
    def __init__(
            self,
            M=10,
            Nstar=4.16,  # [rps] impeller revolutions
            phi=0.15,  # [1] holdup
            vmax=0.08):
        self.D = 10  # [cm] impeller diameter
        self.Nstar = Nstar
        self.phi = phi

        # Water
        self.muc = 0.89e-2  # [P = g * cm^-1 s^-1]
        self.rhoc = 1.0  # [g/cm3]
        # Kerosene-dicholorebenzene
        self.rhod = 0.972  # [g/cm3]
        self.sigma = 42.82  # Unit?

        time = arange(0.0, 3600, 0.5)

        mm3_to_cm3 = 0.1**3  # TODO: Fix that scaling
        vmax = vmax * mm3_to_cm3
        # vmax = 0.06 * mm3_to_cm3

        # Feed distribution
        self.v0 = vmax / 2
        self.sigma0 = (vmax - self.v0) / 3.3
        vmin = None  # 5.23e-7

        # Feed
        theta = 600
        self.Vt = 12 * 10**3
        self.n0 = self.Vt / theta
        MOCSolution.__init__(
            self, M, time, vmax / M, N0=self.N0, xi0=vmin,
            beta=beta, gamma=self.g, Q=self.Qf,
            theta=theta, n0=self.n0, A0=self.A0)

    def N0(self, v):
        return 0 * v

    def A0(self, v):
        return \
            self.phi / (self.v0 * self.sigma0 * sqrt(2 * pi)) * \
            exp(-(v - self.v0)**2 / (2 * self.sigma0**2))

    def g(self, v, C1=0.4, C2=0.08):
        return \
            C1 * v**(-2. / 9) * self.D**(2. / 3) * \
            self.Nstar / (1 + self.phi) * \
            exp(- C2 * (1 + self.phi)**2 * self.sigma /
                (self.rhod * v**(5. / 9) * self.D**(4. / 3) * self.Nstar**2))

    def Qf(self, v1, v2, C3=2.8e-6, C4=1.83e9):
        d_ratio = (v1**(1. / 3) * v2**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))

        return C3 * \
            (v1**(2. / 3) + v2**(2. / 3)) * \
            (v1**(2. / 9) + v2**(2. / 9))**0.5 * \
            self.D**(2. / 3) * \
            self.Nstar / (1 + self.phi) * \
            exp(
                -C4 * self.muc * self.rhoc * self.D**2 /
                self.sigma**2 *
                self.Nstar**3 /
                (1 + self.phi)**3 *
                d_ratio**4)

    @property
    def pbe_phi(self):
        return self.total_volume / self.Vt
