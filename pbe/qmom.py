import numpy as np


def wheeler_inversion(m):
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

    b_diag = -np.sqrt(b[1:])
    jacobi = np.diag(a) + np.diag(b_diag, k=1) + np.diag(b_diag, k=-1)

    eigval, eigvec = np.linalg.eig(jacobi)
    w = m[0] * eigvec[0, :]**2

    return eigval, w

