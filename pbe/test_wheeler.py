import numpy as np
from wheeler import wheeler_inversion
from numpy.testing import assert_almost_equal


#  Friction factor for laminar calculations are 64/Re
#  Applies only for Re lower than 2100
def test_book_ex3_2():
    moments = np.array([1, 0, 1, 0, 3, 0, 15, 0])
    abscissas = np.array([-2.3344, -0.7420, 0.7420, 2.3344])
    weights = np.array([0.0459, 0.4541, 0.4541, 0.0459])
    [xi, w] = wheeler_inversion(moments)
    assert_almost_equal(
        np.max(np.abs(xi-abscissas)), 0,
        decimal=4
    )
    assert_almost_equal(
        np.max(np.abs(w-weights)), 0,
        decimal=4
    )
