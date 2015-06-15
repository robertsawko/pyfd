import fluid
import numpy as np
import matplotlib.pyplot as plt

path = 'validationData/tables/1/'
f = open(path + 'parameters', 'w')
f_br = open(path + 'breakup', 'w')
f_c = open(path + 'coalescence', 'w')

F = fluid.fluid('coulaloglou', 0, 0)
F.C = np.array([1.52e-02, 6.78e-02, 1.06e-12, 5.13e13])

for i in range(4):
    f.write('C{0} '.format(i + 1) + repr(F.C[i]) + '\n')
f.write('rhoc ' + repr(F.rhoc) + '\n')
f.write('rhod ' + repr(F.rhod) + '\n')
f.write('muc ' + repr(F.muc) + '\n')
f.write('epsilon ' + repr(F.epsilon) + '\n')
f.write('sigma ' + repr(F.sigma) + '\n')
f.write('alpha ' + repr(F.alpha) + '\n')
f.write('V ' + repr(F.V) + '\n')
f.close()

dMin = 0.1 * F.expectedD
dMax = 10.0 * F.expectedD
dd = 0.1 * F.expectedD
d = np.arange(dMin, dMax, dd)
v = np.pi * d ** 3 / 6.0
v0 = np.pi * F.expectedD ** 3 / 6.0

gamma = F.gamma(v)
Q = F.Q(v, v0)
v2 = np.zeros(v.shape[0])
v2.fill(v0)

result_breakup = zip(v, gamma)
result_coalescence = zip(v, v2, Q)
np.savetxt(f_br, result_breakup)
np.savetxt(f_c, result_coalescence)

#fig = plt.figure()
#plt.plot(v, gamma)
#ax = fig.gca()
#ax.plot(v, Q)
#ax.set_yscale('log')
#plt.show()
