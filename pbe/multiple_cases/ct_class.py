from case_class import CaseSolution
from numpy import arange, sqrt, exp, pi, zeros


class CTSolution(CaseSolution):
    def __init__(
            self,
            M=10,
            Nstar=4.16,  # [rps] impeller revolutions
            phi=0.15,  # [1] holdup
            v0=0.03):
        self.D = 10  # [cm] impeller diameter

        contProperties = dict()
        dispProperties = dict()
        domainProperties = dict()

        # Water
        contProperties['mu'] = 0.89e-2  # [P = g * cm^-1 s^-1]
        contProperties['rho'] = 0.1  # [g/cm3]
        contProperties['epsilon'] = 0.407 * Nstar * self.D
        # Kerosene-dicholorebenzene
        dispProperties['sigma'] = 42.82  # [P = g * cm^-1 s^-1]
        dispProperties['rho'] = 0.972  # [g/cm3]
        dispProperties['phi'] = phi


        mm3_to_cm3 = 0.1**3
        dispProperties['vMax'] = 0.08 * mm3_to_cm3

        # Feed distribution
        dispProperties['v0'] = v0 * mm3_to_cm3
        dispProperties['sigma0'] = 0.005 * mm3_to_cm3

        # Feed
        domainProperties['theta'] = 600.
        domainProperties['V'] = 12 * 10**3
        domainProperties['M'] = M
        Ninit = zeros(M)
        CaseSolution.__init__(
            self, dispProperties, contProperties, domainProperties)
