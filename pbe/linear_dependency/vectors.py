import numpy as np

# variables
# epsilon, rhoc, rhod, d32, sigma, U, L, muc
Re = np.array([0, 1, 0, 0, 0, 1, 1, -1])
Ca = np.array([0, 0, 0, 0, -1, 1, 0, 1])
St = np.array([0, 0, 0, 2, 0, 0, -2, 0]) + Re
We = np.array([2./3., 1, 0, 5./3., -1, 0, 0, 0])

A = np.array([Re, St, Ca, We])

# if Re, Ca, St and We are linearly independent this should read 4
print np.linalg.matrix_rank(A)
