import numpy as np
from scipy.optimize import fsolve


def reynolds_number(U, L, nu):
    return U * L / nu


def colebrook_white_equation(f, Re, Dh, epsilon=0):
    """
    colebrook_white_equation description

    @type f: float
    @param f: friction factor value

    @type epsilon: float
    @param epsilon: dimensional rougness parameter

    @type Re: float
    @param Re: Reynolds number

    @type Dh: float
    @param Dh: Hydraulic diameter
    """
    return 1 / np.sqrt(f) + \
        2 * np.log10(2.51 / (Re * np.sqrt(f) + epsilon / (3.7 * Dh)))


def friction_factor(Re, Dh):
    if Re > 2100:
        def friction_eq(f):
            return colebrook_white_equation(f, Re, Dh)
        f = fsolve(friction_eq, 64. / Re)
        return f
    else:
        return 64. / Re


def pressure_drop(U, L, D, rho, nu):
    Re = reynolds_number(U, D, nu)
    f = friction_factor(Re, D)
    return rho * f * L / D * (U**2 / 2)


# Geometry specification
pipe_diameter = 48 * 0.0254
orifice_diameter = pipe_diameter / 12.5

# Oil properties
oil_kinematic_viscosity = 1.9484e-05
oil_density = 871.5222

# Flow parameters
oil_inlet_velocity = 0.3  # m/s
pig_velocity = 0.2  # m/s

orifice_velocity = (
    (oil_inlet_velocity - pig_velocity)
    * (pipe_diameter / orifice_diameter)**2)

print orifice_velocity

dp = pressure_drop(
    orifice_velocity, pipe_diameter, orifice_diameter,
    oil_density, oil_kinematic_viscosity)

print dp
