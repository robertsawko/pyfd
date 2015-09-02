from numpy import arange, linspace, zeros, array, sqrt
from itertools import cycle
from moc import MOCSolution
from test_moc import ziff_total_number_solution, ziff_pbe_solution
import matplotlib as mpl
mpl.use('pgf')
import matplotlib.pyplot as plt
plt.style.use('ggplot')


def set_plt_params(
        relative_fig_width=1.0, landscape=True, page_width=307.3, rescale_h=1):

    fig_width_pt = page_width * relative_fig_width
    inches_per_pt = 1.0 / 72.27               # Convert pt to inch
    golden_mean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
    fig_width = fig_width_pt * inches_per_pt  # width in inches

    if landscape:
        fig_height = fig_width * golden_mean  # height in inches
    else:
        fig_height = fig_width / golden_mean  # height in inches

    fig_height = fig_height * rescale_h
    fig_size = [fig_width, fig_height]
    params = {
        'font.family': 'serif',
        'axes.labelsize': 7,
        'xtick.labelsize': 5,
        'ytick.labelsize': 5,
        'axes.labelcolor': 'black',
        'ytick.color': 'black',
        'xtick.color': 'black',
        'legend.handlelength': 4,
        'legend.fontsize': 7,
        # 'lines.markersize': 3,
        # 'xtick.labelsize': 7,
        # 'ytick.labelsize': 7,
        'text.usetex': True,
        'text.latex.unicode': True,
        'figure.figsize': fig_size,
        'pgf.texsystem': "xelatex",
        'pgf.rcfonts': False,
    }

    plt.rcParams.update(params)
"""
Case setup based on:

    Ziff, R.M and McGrady, E.D
    "New solutions to the fragmentation equation", J. Phys A, vol. 18, 1985

    We are looking at case 4. from their paper with kernel F(x, y) = x + y.
    This corresponds to choosing:
        beta = 2.0/y,
        Gamma = y^2,
    for our kernels.
"""

grids = [20, 40, 80, 160]
time = arange(0.0, 10.0, 0.001)
vmax = 1.0

pbe_solutions = dict()
for g in grids:
    N0 = zeros(g)
    N0[-1] = 1
    pbe_solutions[g] = MOCSolution(
        N0, time, vmax / g,
        beta=lambda x, y: 2.0 / y,
        gamma=lambda x: x**2
    )

totals = dict(
    (
        n,
        array([sum(Ns) for Ns in pbe_solutions[n].N])
    ) for n in pbe_solutions
)

v = linspace(0, vmax, 100)
Na = array([ziff_total_number_solution(v, t, vmax) for t in time])


set_plt_params(relative_fig_width=0.49)
fig = plt.figure()
ax = fig.gca()
linestyles = cycle(['-', '--', '-.', ':'])
for n in sorted(totals):
    ax.loglog(
        time, totals[n] / totals[n][0],
        linestyle=next(linestyles),
        label='MOC with N={0}'.format(n))
ax.loglog(time, Na, "-k", linewidth=2, label='Analytical')
ax.set_xlabel('Time $t$')
ax.set_ylabel('Total number of drops $N/N_0$')
plt.savefig("pure_breakup-total_number.pgf", tight_layout=True)


set_plt_params(relative_fig_width=0.49)
fig = plt.figure()
ax = fig.gca()
markers = cycle(['o', 's', 'v', '*', '.', ','])
v = linspace(0, vmax, 10000)

for n in sorted(pbe_solutions):
    ax.semilogy(
        pbe_solutions[n].xi, pbe_solutions[n].number_density()[-1], "+",
        marker=next(markers),
        label="MOC with N={0}".format(n))
ax.semilogy(
    v, ziff_pbe_solution(v, time[-1], vmax), "-k",
    linewidth=2, label="Analytical $t=\infty$")
ax.legend(loc='lower left', shadow=True)
ax.set_xlabel('Particle volume')
ax.set_ylabel('Number density function')
plt.savefig("pure_breakup-pbe.pgf", tight_layout=True)
figlegend = ax.legend(shadow=True)
