import numpy as np
from scipy.optimize import fsolve


def reynolds_number(U, L, nu):
    return U * L / nu


def colebrook_white_equation(f, Re, Dh, epsilon=0):
    return (
        1/np.sqrt(f)
        + 2 * np.log10(2.51 / (Re*np.sqrt(10) + epsilon/(3.7 * Dh))))


def friction_factor(Re, Dh):
    def friction_eq(f):
        return colebrook_white_equation(f, Re, Dh)

    f = fsolve(friction_eq, 64 / Re)
    return f

# Geometry specification
pipe_diameter = 48 * 0.0254
orifice_diameter = pipe_diameter / 12.5

# Oil properties
oil_kinematic_viscosity = 1.9484e-05
oil_density = 871.5222

# Flow parameters
oil_inlet_velocity = 0.3  # m/s
pig_velocity = 0.0  # m/s

Re_pipe = reynolds_number(
    oil_inlet_velocity, pipe_diameter, oil_kinematic_viscosity)

orifice_velocity = (
    (oil_inlet_velocity - pig_velocity)
    * (pipe_diameter / orifice_diameter)**2)

Re_orifice = reynolds_number(
    orifice_velocity, orifice_diameter, oil_kinematic_viscosity)

print Re_pipe
print Re_orifice
print friction_factor(Re_pipe, pipe_diameter)
