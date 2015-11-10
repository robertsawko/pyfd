from numpy import genfromtxt, abs, array, pi
from scipy.optimize import minimize
from pyfd.pbe.moc import CTSolution
import time
import pickle


class ct_experiment:
    def __init__(self, phi, Nstar, d32):
        self.phi = phi
        self.Nstar = Nstar
        self.d32 = d32

def error_function(C, experiment):
    v0s = array([0.5, 1.5]) * pi / 6 * experiment.d32**3
    mp = array(C)
    mp[3] *= 1e13
    pbe_solutions = [
        CTSolution(
            M=20, v0=v0, Nstar=experiment.Nstar, phi=experiment.phi,
            model_parameters=mp)
        for v0 in v0s]
    error = sum(
        [abs(s.d32 - experiment.d32) / experiment.d32
            for s in pbe_solutions])
    print('constants={0}, error={1}'.format(C, error))
    return error

concentrations = [5, 10, 15]
mm_to_m = 0.001

experiments = []
for c in concentrations:
    data = genfromtxt(
        '../../../pyfd/data/pbe/coulaloglou/d32_N_alpha{0}.txt'.format(c))
    for i in range(data.shape[0]):
        experiments.append(
            ct_experiment(c / 100., data[i, 0] / 60, data[i, 1] * mm_to_m))


c0 = [0.4, 0.08, 2.8, 1.83]  # CT original constants
results = []
for e in experiments[0:1]:
    res = dict()

    Copt = minimize(
        lambda c: error_function(c, e), c0,
        method='TNC', bounds=[(0, None)] * 4,
        options={'disp': False, 'ftol': 0.01, 'maxiter': 1})
    res['best_fit'] = Copt.x
    res['setup'] = {'phi': e.phi, 'Nstar': e.Nstar, 'd32': e.d32}
    results.append(res)

with open('data/ct_optimization_results.pickle', 'wb') as f:
    pickle.dump(results, f)
