import numpy as np
from pyfd.pbe.qmom import is_realizable
from numpy.testing import assert_equal

'''
Unit tests for realizability

Based on example 3.3

Daniele L. Marchisio, Rodney O. Fox "Computational Models for Polydisperse
Particulate and Multiphase Systems", Cambridge University Press 2013
'''


def test_book_ex3_3():
    moments = np.array([1, 5, 26, 140, 778, 4450, 26140, 157400])
    assert_equal(is_realizable(moments), True)
    moments[1] = 26
    assert_equal(is_realizable(moments), False)
