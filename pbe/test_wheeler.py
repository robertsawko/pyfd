import numpy as np
from wheeler import wheeler_inversion
from numpy.testing import assert_almost_equal


#  Friction factor for laminar calculations are 64/Re
#  Applies only for Re lower than 2100
def test_book_ex3_2():
    moments = np.array([1, 0, 1, 0, 3, 0, 15, 0])
    xi1 = np.sqrt(np.sqrt(6) + 3)  # Taken from Maxima
    xi2 = np.sqrt(3 - np.sqrt(6))
    w1 = 1/(2*np.sqrt(6) + 2**(3./2)*np.sqrt(3) + 12)
    w2 = 1/(12 - 2*np.sqrt(6) - 2**(3./2)*np.sqrt(3))
    abscissas = np.array([-xi1, -xi2, xi2, xi1])
    weights = np.array([w1, w2, w2, w1])
    [xi, w] = wheeler_inversion(moments)
    assert_almost_equal(np.max(np.abs(xi-abscissas)), 0)
    assert_almost_equal(np.max(np.abs(w-weights)), 0)
