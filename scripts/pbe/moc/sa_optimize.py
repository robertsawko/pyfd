from numpy import genfromtxt, abs, array, pi
from scipy.optimize import minimize, differential_evolution
from pyfd.pbe.moc import SASolution
import time
import pickle


class sa_experiment:
    def __init__(self, d32=0.32e-03, phi=0.117, U=2.72, theta=600.):
        self.phi = phi
        self.theta = theta
        self.U = U
        self.d32 = d32

def error_function(C, experiment):
    v0s = array([0.5, 1.5]) * pi / 6 * experiment.d32**3
    mp = array(C)
    mp[3] *= 1e13
    pbe_solutions = [
        SASolution(
            M=20, v0=v0, U=experiment.U, phi=experiment.phi, theta=experiment.theta,
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
experiments.append(sa_experiment())
# original C&T parameters where multiplying (Nstar**3 * D**2)
# and epsilon = 0.407 * NStar**3 * D**2
# so
s=0.407
c0 = [0.4 * s**(-1./3.), 0.08 / s**(-2./3.), 2.8 * s**(-1./3.), 1.83 * s]  # CT original constants
#c0 = [0.75, 0.22, 5., 4.3]
results = []
for e in experiments:
    res = dict()

    Copt = minimize(
        lambda c: error_function(c, e), c0,
        method='L-BFGS-B', bounds=[(0, None)] * 4,
        options={'disp': False, 'ftol': 0.01, 'maxiter': 100})

    #b = [(0, 1), (0, 1), (0, 10), (0, 10)] 
    #Copt = differential_evolution(lambda c: error_function(c, e), 
        #bounds=b, maxiter=10, popsize=1)
    res['best_fit'] = Copt.x
    res['setup'] = {'U': e.U, 'phi': e.phi, 'd32': e.d32, 'theta': e.theta}
    results.append(res)

with open('data/sa_optimization_results.pickle', 'wb') as f:
    pickle.dump(results, f)
