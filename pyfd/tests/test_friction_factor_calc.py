from pyfd.calc.friction_factor import friction_factor, reynolds_number
from pyfd.calc.friction_factor import pressure_drop
from numpy.testing import assert_approx_equal


#  Friction factor for laminar calculations are 64/Re
#  Applies only for Re lower than 2100
def test_laminar():
    for Re in [1, 10, 100, 1000]:
        assert_approx_equal(
            friction_factor(Re, 1),
            64./Re,
            significant=2
        )


D = 48 * 0.0254
U = 0.3
nu = 1.9484e-05
rho = 871.5222
L = 50


#  Numerical results based on
#  http://www.efunda.com/formulae/fluids/calc_pipe_friction.cfm#calc
def test_f_vs_efunda():
    assert_approx_equal(
        friction_factor(reynolds_number(U, D, nu), D),
        0.0263,
        significant=2
    )


def test_dp_vs_efunda():
    assert_approx_equal(
        pressure_drop(U, L, D, rho, nu),
        0.0422 * 1000,  # 0.0422 kPa
        significant=2
    )
