from case_class import CaseSolution, DomainProperties, DispersionProperties
from numpy import pi


class AngeliSolution(CaseSolution):
    def __init__(
            self,
            M=10,
            U=1.1,  # [rps] impeller revolutions
            phi=0.091,  # [1] holdup
            v0=5e-10,  # [cm^3]
            model_parameters=None,
            theta=600.):

        # pipe diamter and length
        self.D = 24.0e-03
        self.L = 9.5

        contProperties = dict()

        # oil
        contProperties['mu'] = 1.6e-3  # [P = kg * m^-1 s^-1]
        contProperties['rho'] = 801.  # [kg/cm3]

        # calculate turbulent properties
        Re = U * self.D / contProperties['mu'] * contProperties['rho']
        I = 0.16 * Re ** (-1. / 8.)
        u_rms = U * I
        k = 3. / 2. * u_rms ** 2
        L_t = 0.038 * self.D
        contProperties['epsilon'] = 0.09 * k ** (3. / 2.) / L_t
        contProperties['Re'] = Re
        # water solution
        dispersion = DispersionProperties(
            phi=phi, rho=1000.,  # [kg/m3]
            sigma=1.7e-2,  # [P = kg * m^-1 s^-1]
            v_max=v0 * 3,
            v0=v0,
            sigma0=v0 / 10)

        # Feed distribution

        # Feed
        domain = DomainProperties(
            theta=theta, V=pi * self.L * (self.D / 2) ** 2, M=M)

        CaseSolution.__init__(
            self, dispersion, contProperties, domain,
            model_parameters=model_parameters)
