from friction_factor_calc import pressure_drop
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
    (oil_inlet_velocity - pig_velocity) *
    (pipe_diameter / orifice_diameter)**2)

print(orifice_velocity)

dp = pressure_drop(
    orifice_velocity, pipe_diameter, orifice_diameter,
    oil_density, oil_kinematic_viscosity)

print(dp)
