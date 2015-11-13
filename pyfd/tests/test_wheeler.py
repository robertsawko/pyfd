import numpy as np
from pyfd.pbe.qmom import wheeler_inversion
from numpy.testing import assert_almost_equal

"""
Unit tests for Wheeler moment inversion algorithm.

Based on examples 3.1 and 3.2 from:

Daniele L. Marchisio, Rodney O. Fox "Computational Models for Polydisperse
Particulate and Multiphase Systems", Cambridge University Press 2013

The book only gives the results up to for decimal digit precision. The
analytical results were obtained using Maxima by solving the Jacobians 3.18 and
3.24. The second case is interesting as the product difference algorithm
wouldn't work there.
"""


def test_book_ex3_1():
    moments = np.array([1, 5, 26, 140, 778, 4450, 26140, 157400])
    w1 = 1/(2*np.sqrt(6) + 2**(3./2)*np.sqrt(3) + 12)
    w2 = 1/(12 - 2*np.sqrt(6) - 2**(3./2)*np.sqrt(3))
    abscissas = np.array([2.6656, 4.2580, 5.7420, 7.3344])
    weights = np.array([w1, w2, w2, w1])
    [xi, w] = wheeler_inversion(moments)
    assert_almost_equal(np.max(np.abs(xi-abscissas)), 0, decimal=4)
    assert_almost_equal(np.max(np.abs(w-weights)), 0)


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
