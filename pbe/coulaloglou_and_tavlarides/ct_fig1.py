from ct_class import CTSolution
import matplotlib.pyplot as plt

sol = CTSolution(M=750, Nstar=310. / 60, phi=0.1)
plt.plot(sol.xi_d * 10, sol.number_density[2800])
plt.plot(sol.xi_d * 10, sol.number_density[-1])
plt.show()
