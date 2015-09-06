from ct_class import CTSolution
from numpy import linspace
import pickle

concentrations = [5, 10, 15]
Ns = linspace(3, 6, 10)  # rps
ct_solutions = dict([(
    c, [CTSolution(M=20, Nstar=N, phi=c / 100.0) for N in Ns])
    for c in concentrations])

with open('ct_solutions.pickle', 'wb') as f:
    pickle.dump(ct_solutions, f)
