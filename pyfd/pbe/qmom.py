"""
Functions for quadrature based method of moments

0. Orthogonal polynomial theory
1. Moment inversion algorithms
2. Moment realizability checks
"""
import numpy as np


def calculate_weights_and_abcissas(a, b, m_0):
    b_diag = -np.sqrt(np.abs(b))
    jacobi = np.diag(a) + np.diag(b_diag, k=1) + np.diag(b_diag, k=-1)

    eigval, eigvec = np.linalg.eig(jacobi)
    w = m_0 * eigvec[0, :]**2

    return eigval, w


def product_difference_inversion(m):
    """Return abcissas and weights for a given set of moments using product
    difference algorithm.

    :param m: array of moments of length N
    :returns: weights and abcissas :math:`(\{xi_i\}_{i=0}^N,\{w_i\}_{i=0}^N)`
    :rtype: Tuple(array, array)
    """
    N = int(m.size / 2)

    Pn = 2 * N + 1

    P = np.zeros((Pn, Pn))
    P[0, 0] = 1
    for i in range(Pn-1):
        P[i, 1] = (-1)**i*m[i]
    for j in range(2, Pn):
        for i in range(Pn+1-j):
            P[i, j] = P[0, j-1]*P[i+1, j-2]-P[0, j-2]*P[i+1, j-1]

    zeta = np.concatenate(([0], P[0, 2:2*N+1]/(P[0, 1:2*N]*P[0, :2*N-1])))

    a = zeta[0::2] + zeta[1::2]
    b = zeta[1:-1:2] * zeta[2::2]

    return calculate_weights_and_abcissas(a, b, m[0])


def wheeler_inversion(m):
    """Return abcissas and weights for a given set of moments using Wheeler
    algorithm (Chebyshev algorithm).

    :param m: array of moments of length N
    :returns: weights and abcissas :math:`(\{xi_i\}_{i=0}^N,\{w_i\}_{i=0}^N)`
    :rtype: Tuple(array, array)
    """
    N = int(m.size / 2)  # number of nodes of the quadrature approximation

    sigma = np.zeros((N+1, 2*N))
    sigma[1, :] = m[0:2*N]

    a = np.zeros(N)
    b = np.zeros(N)
    a[0] = m[1] / m[0]
    b[0] = m[0]  # This value is insignificant as it's not being used.

    for k in range(1, N):
        l = np.arange(2*(N - k)) + k
        sigma[k+1, l] = sigma[k, l+1] - a[k-1]*sigma[k, l] - \
            b[k-1]*sigma[k-1, l]

        a[k] = -sigma[k, k] / sigma[k, k-1] + \
            sigma[k+1, k+1] / sigma[k+1, k]

        b[k] = sigma[k+1, k] / sigma[k, k-1]

    return calculate_weights_and_abcissas(a, b[1:], m[0])


def is_realizable(m):
    N = int(m.size/2)
    M = np.zeros((N, N))
    for i in range(0, N):
        M[i, :] = m[i:N+i]
    for l in range(1, N):
        if (np.linalg.det(M[0:l, 0:l]) < 0):
            return False
    for l in range(2, N-1):
        if (np.linalg.det(M[1:l, 1:l]) < 0):
            return False

    return True
