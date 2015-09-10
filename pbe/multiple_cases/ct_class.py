from case_class import CaseSolution
from numpy import zeros


class CTSolution(CaseSolution):
    def __init__(
            self,
            M=10,
            Nstar=4.16,  # [rps] impeller revolutions
            phi=0.15,  # [1] holdup
            v0=3e-11,  # [cm^3]
            model_parameters=None):
        self.D = 0.10  # [m] impeller diameter

        contProperties = dict()
        dispProperties = dict()
        domainProperties = dict()

        # Water
        contProperties['mu'] = 0.89e-3  # [P = kg * m^-1 s^-1]
        contProperties['rho'] = 1000  # [kg/cm3]
        # contProperties['epsilon'] = 0.407 * Nstar**3 * self.D**2
        contProperties['epsilon'] = Nstar**3 * self.D**2
        # Kerosene-dicholorebenzene
        dispProperties['sigma'] = 42.82e-3  # [P = kg * m^-1 s^-1]
        dispProperties['rho'] = 972  # [kg/m3]
        dispProperties['phi'] = phi

        dispProperties['vMax'] = 6e-11

        # Feed distribution
        dispProperties['v0'] = v0
        dispProperties['sigma0'] = v0 / 10

        # Feed
        domainProperties['theta'] = 600.
        domainProperties['V'] = 12e-3
        domainProperties['M'] = M
        CaseSolution.__init__(
            self, dispProperties, contProperties, domainProperties,
            model_parameters=model_parameters)
