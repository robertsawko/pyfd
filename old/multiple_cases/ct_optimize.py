from numpy import genfromtxt, abs, array, pi
from scipy.optimize import minimize
from ct_class import CTSolution


class ct_experiment:
    def __init__(self, phi, Nstar, d32):
        self.phi = phi
        self.Nstar = Nstar
        self.d32 = d32

concentrations = [5, 10, 15]

mm_to_cm = 0.1
experiments = []
for c in concentrations:
    data = genfromtxt(
        '../validationData/coulaloglou/d32_N_alpha{0}.txt'.format(c))
    for i in range(data.shape[0]):
        experiments.append(
            ct_experiment(c / 100., data[i, 0] / 60, data[i, 1] * mm_to_cm))


def error_function(C, experiment):
    v0s = array([0.8, 1.2]) * pi / 6 * experiment.d32**3
    pbe_solutions = [
        CTSolution(
            M=10, v0=v0, Nstar=experiment.Nstar, phi=experiment.phi,
            model_parameters=C)
        for v0 in v0s]
    error = sum(
        [abs(s.d32 - experiment.d32) / experiment.d32
            for s in pbe_solutions])
    print('constants={0}, error={1}'.format(C, error))
    return error

c0 = [0.4, 0.08, 2.8e-6, 1.83e9]  # CT original constants
best_fit_constants = []
for e in experiments[0:1]:
    # do our magic optimization
    #print(
        #'Optimizing phi={0:0.2f}, Nstar={1:0.2f}'.format(e.phi, e.Nstar),
        #end='...', flush=True)
    Copt = minimize(
        lambda c: error_function(c, e), c0,
        method='TNC', bounds=[(0, None)] * 4,
        options={'disp': False})
    #if Copt.success:
        #print('success.')
    #else:
        #print('failed.')
    best_fit_constants.append(Copt.x)
