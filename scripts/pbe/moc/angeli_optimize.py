from numpy import genfromtxt, abs, array, pi
from scipy.optimize import minimize, differential_evolution
from pyfd.pbe.moc import AngeliSolution
import time
import pickle


class angeli_experiment:
    def __init__(self, U, phi, d32, theta=1200.):
        self.phi = phi
        self.theta = theta
        self.U = U
        self.d32 = d32

def error_function(C, experiment):
    v0s = array([0.5, 1.5]) * pi / 6 * experiment.d32**3
    mp = array(C)
    mp[3] *= 1e11
    pbe_solutions = [
        AngeliSolution(
            M=20, v0=v0, U=experiment.U, phi=experiment.phi, theta=experiment.theta,
            model_parameters=mp)
        for v0 in v0s]
    error = sum(
        [abs(s.d32 - experiment.d32) / experiment.d32
            for s in pbe_solutions])
    print('constants={0}, error={1}'.format(C, error))
    print('d=[{0}, {1}], expd={2}'.format(pbe_solutions[0].d32, pbe_solutions[1].d32, experiment.d32))
    return error

concentrations = [5, 10, 15]
mm_to_m = 0.001

experiments = []
Us = [1.1, 1.25, 1.3, 1.5, 1.7]
alphas = [0.091, 0.05, 0.077, 0.034, 0.05]
ds = [1.29, 0.97, 1.25, 0.73, 0.69]
for i in range(5):
    experiments.append(angeli_experiment(Us[i], alphas[i], ds[i] * 1e-03))

# original C&T parameters where multiplying (Nstar**3 * D**2)
# and epsilon = 0.407 * NStar**3 * D**2
# so
s=0.407
c0 = [0.4 * s**(-1./3.), 0.08 / s**(-2./3.), 2.8 * s**(-1./3.), 1.83 * s]  # CT original constants

results = []
for e in experiments:
    res = dict()

    Copt = minimize(
        lambda c: error_function(c, e), c0,
        method='L-BFGS-B', bounds=[(0, None)] * 4,
        options={'disp': False, 'ftol': 0.01, 'maxiter': 100})
    res['best_fit'] = Copt.x
    res['setup'] = {'U': e.U, 'phi': e.phi, 'd32': e.d32, 'theta': e.theta}
    results.append(res)

with open('data/angeli_optimization_results.pickle', 'wb') as f:
    pickle.dump(results, f)
