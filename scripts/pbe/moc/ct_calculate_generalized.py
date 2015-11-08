from pyfd.pbe.moc import CTSolution
from numpy import linspace
import pickle

"""
This script calculated solutions to PBE problem that are necessary to reproduce
figure 3 from CT publication.
"""
concentrations = [5]  # , 10, 15]
Ns = linspace(3, 6, 10)  # rps
ct_solutions = dict([(
    c, [CTSolution(M=50, Nstar=N, phi=c / 100.0) for N in Ns])
    for c in concentrations])

with open('data/ct_solutions_generalized.pickle', 'wb') as f:
    pickle.dump(ct_solutions, f)
