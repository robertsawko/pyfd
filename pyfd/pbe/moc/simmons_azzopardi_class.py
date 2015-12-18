from case_class import CaseSolution, DomainProperties, DispersionProperties
from numpy import pi


class SASolution(CaseSolution):
    def __init__(
            self,
            M=10,
            U=2.71,  # [rps] impeller revolutions
            phi=0.117,  # [1] holdup
            v0=4e-11,  # [cm^3]
            model_parameters=None,
            theta=600.):

        # pipe diamter andlength
        self.D = 0.063  # [m] impeller diameter
        self.L = 4.5    # [m] impeller diameter

        contProperties = dict()

        # oil
        contProperties['mu'] = 1.8e-3  # [P = kg * m^-1 s^-1]
        contProperties['rho'] = 797.  # [kg/cm3]

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
            phi=phi, rho=1166.,  # [kg/m3]
            sigma=1.e-2,  # [P = kg * m^-1 s^-1]
            v_max=6e-11,
            v0=v0,
            sigma0=v0 / 10)

        domain = DomainProperties(
            theta=theta, V=pi * self.L * (self.D / 2) ** 2, M=M)

        CaseSolution.__init__(
            self, dispersion, contProperties, domain,
            model_parameters=model_parameters)
