import tables
from itertools import izip


def plot_ldf_and_models(data, group):
    global binned_densities
    clf()
    particles = data.get_node(group)

    bins = logspace(1, 3, 50)
    x = (bins[:-1] + bins[1:]) / 2

    # plot ldf from ground particles
    R = particles[:]['core_distance']
    N, hbins = histogram(R, bins)
    area = [pi * (v ** 2 - u ** 2) for u, v in zip(bins[:-1], bins[1:])]
    loglog(x, N / area, label="full")

    # plot ldf from observables
    observables = data.root.simulations.E_1PeV.zenith_0.observables
    R = observables[:]['r']
    n1 = observables[:]['n1']
    n2 = observables[:]['n2']
    n3 = observables[:]['n3']
    n4 = observables[:]['n4']
    density = (n1 + n2 + n3 + n4) / 2.

    R_within_limits = R.compress((R < bins[-1]) & (R >= bins[0]))
    density_within_R_limits = density.compress((R < bins[-1]) & (R >= bins[0]))

    idxs = searchsorted(bins, R_within_limits) - 1
    binned_densities = [[] for i in range(len(bins) - 1)]
    for idx, value in zip(idxs, density_within_R_limits):
        binned_densities[idx].append(value)
    y = [mean(u) for u in binned_densities]
    y_err = [std(u) for u in binned_densities]
    errorbar(x, y, y_err, label="measured")

    xlabel("Core distance [m]")
    ylabel("Particle density [m^{-2}]")
    title("Full and measured LDF (E = 1 PeV)")
    legend()


def plot_ldf_ldf(data, group):
    group = data.get_node(group)

    bins = logspace(1, 3, 50)
    x = (bins[:-1] + bins[1:]) / 2

    # plot ldf from observables
    R, density = [], []
    for c in group.coincidences[:10000]:
        shower_x = c['x']
        shower_y = c['y']
        for obs_idx in group.c_index[c['id']]:
            observables = group.observables[obs_idx]
            obs_x = observables['x']
            obs_y = observables['y']
            R.append(sqrt((obs_x - shower_x) ** 2 + (obs_y - shower_y) ** 2))
            n1 = observables['n1']
            n2 = observables['n2']
            n3 = observables['n3']
            n4 = observables['n4']
            density.append((n1 + n2 + n3 + n4) / 2.)

    R = array(R)
    density = array(density)

    R_within_limits = R.compress((R < bins[-1]) & (R >= bins[0]))
    density_within_R_limits = density.compress((R < bins[-1]) & (R >= bins[0]))

    idxs = searchsorted(bins, R_within_limits) - 1
    binned_densities = [[] for i in range(len(bins) - 1)]
    for idx, value in zip(idxs, density_within_R_limits):
        binned_densities[idx].append(value)
    y = [mean(u) for u in binned_densities]
    y_err = [std(u) for u in binned_densities]
    errorbar(x, y, y_err, label="measured ldf_sim")
    legend()


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.open_file('data-e15-S250.h5', 'r')

    sim = data.root.simulations.E_1PeV.zenith_0
    plot_ldf_and_models(data, '/showers/E_1PeV/zenith_0/electrons')
    plot_ldf_ldf(data, '/ldfsim')
    savefig('plots/ldfs_sim_1PeV.pdf')
