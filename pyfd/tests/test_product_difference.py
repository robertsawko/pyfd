import numpy as np
from pyfd.pbe.qmom import product_difference_inversion
from numpy.testing import assert_almost_equal


def test_book_ex3_1():
    moments = np.array([1, 5, 26, 140, 778, 4450, 26140, 157400])
    w1 = 1/(2*np.sqrt(6) + 2**(3./2)*np.sqrt(3) + 12)
    w2 = 1/(12 - 2*np.sqrt(6) - 2**(3./2)*np.sqrt(3))
    abscissas = np.array([2.6656, 4.2580, 5.7420, 7.3344])
    weights = np.array([w1, w2, w2, w1])
    [xi, w] = product_difference_inversion(moments)
    assert_almost_equal(np.max(np.abs(xi-abscissas)), 0, decimal=4)
    assert_almost_equal(np.max(np.abs(w-weights)), 0)
