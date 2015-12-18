from .case_class import CaseSolution, DomainProperties, DispersionProperties


class CTSolution(CaseSolution):
    def __init__(
            self,
            M=10,
            Nstar=4.16,  # [rps] impeller revolutions
            phi=0.15,  # [1] holdup
            v0=4e-11,  # [cm^3]
            model_parameters=None):
        self.D = 0.10  # [m] impeller diameter

        contProperties = dict()

        # Water
        contProperties['mu'] = 0.89e-3  # [P = kg * m^-1 s^-1]
        contProperties['rho'] = 1000  # [kg/cm3]
        # contProperties['epsilon'] = 0.407 * Nstar**3 * self.D**2
        contProperties['epsilon'] = Nstar**3 * self.D**2
        # Kerosene-dicholorebenzene

        dispersion = DispersionProperties(
            phi=phi, rho=972.,  # [kg/m3]
            sigma=42.82e-3,  # [P = kg * m^-1 s^-1]
            v_max=6e-11,
            v0=v0,
            sigma0=v0 / 10.)

        domain = DomainProperties(theta=600, V=12e-3, M=M)

        CaseSolution.__init__(
            self, dispersion, contProperties, domain,
            model_parameters=model_parameters)
