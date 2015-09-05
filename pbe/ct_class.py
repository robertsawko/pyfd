from moc import MOCSolution
from numpy import arange, sqrt, exp, pi, zeros


def beta(v1, v2):
    return 2. * 2.4 / v2 * exp(-4.5 * (2 * v1 - v2)**2 / (v2**2))


class CTSolution(MOCSolution):
    def __init__(self, M=10, v0=0.01, n0=1.):
        self.D = 10  # [cm] impeller diameter
        self.Nstar = 4.16  # [rps] impeller revolutions
        self.phi = 0.15  # [1] holdup

        # Water
        self.muc = 1.0e-2  # [P = g * cm^-1 s^-1]
        self.rhoc = 0.1  # [g/cm3]
        # Kerosene-dicholorebenzene
        self.rhod = 0.972  # [g/cm3]
        self.sigma = 42.82

        time = arange(0.0, 3600, 1)
        # time = arange(0.0, 0.02, 0.01)

        mm3_to_cm3 = 0.1**3
        vmax = 0.2 * mm3_to_cm3

        # Feed distribution
        self.v0 = v0 * mm3_to_cm3
        sigma0 = 0.005 * mm3_to_cm3

        # Feed
        theta = 600
        self.Vt = 12 * 10**3
        self.n0 = n0 * self.Vt / theta
        Ninit = zeros(M)
        Ninit[-1] = 1
        dv = vmax / M
        v = dv + dv * arange(M)
        A0 = self.phi / (self.v0 * sigma0 * sqrt(2 * pi)) * \
            exp(-(v - self.v0)**2 / (2 * sigma0**2)) * dv
        MOCSolution.__init__(
            self, Ninit, time, dv,
            beta=beta, gamma=self.g, Q=self.Q,
            theta=theta, n0=self.n0, A0=A0)

    def g(self, v, C1=0.4, C2=0.08):
        return \
            C1 * v**(-2. / 9) * self.D**(2. / 3) * \
            self.Nstar / (1 + self.phi) * \
            exp(- C2 * (1 + self.phi)**2 * self.sigma /
                (self.rhod * v**(5. / 9) * self.D**(4. / 3) * self.Nstar**2))

    def Q(self, v1, v2, C3=2.8e-6, C4=1.83e9):
        ratio = (v1**(1. / 3) * v2**(1. / 3)) / (v1**(1. / 3) + v2**(1. / 3))

        return C3 * \
            (v1**(2. / 3) + v2**(2. / 3)) * \
            (v1**(2. / 9) + v2**(2. / 9))**0.5 * \
            self.D**(2. / 3) * \
            self.Nstar / (1 + self.phi) * \
            exp(
                -C4 * self.muc * self.rhoc * self.D**2 / self.sigma**2 *
                self.Nstar**3 / (1 + self.phi)**3 * ratio**4)

    @property
    def pbe_phi(self):
        return self.total_volume / self.Vt

