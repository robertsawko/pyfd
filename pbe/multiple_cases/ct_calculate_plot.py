from numpy import genfromtxt, linspace
import pickle
import matplotlib.pyplot as plt

with open('ct_solutions_mc.pickle', 'rb') as f:
    ct_solutions = pickle.load(f)

concentrations = [5, 10, 15]
Ns = linspace(3, 6, 10)  # rps

ct_data = dict([(
    c, genfromtxt(
        '../validationData/coulaloglou/d32_N_alpha{0}.txt'.format(c)))
    for c in concentrations])

for c in concentrations:
    plt.plot(ct_data[c][:, 0], ct_data[c][:, 1], 's', )

for c in concentrations:
    plt.plot(Ns * 60, [cts.d32*10 for cts in ct_solutions[c]])

plt.show()
