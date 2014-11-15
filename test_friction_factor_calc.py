from friction_factor_calc import friction_factor, reynolds_number
from nose.tools import assert_almost_equal


#  Results based on
#  http://www.efunda.com/formulae/fluids/calc_pipe_friction.cfm#calc
def test_vs_efunda():
    D = 48 * 0.0254
    U = 0.3
    nu = 1.9484e-05
    print U*D/nu
    assert_almost_equal(
        friction_factor(reynolds_number(U, D, nu), D),
        0.0263,
        delta=1e-4
    )
