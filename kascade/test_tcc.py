import tables
import time

import progressbar as pb

import artist


X, Y = 65., 20.82


progressbar = lambda x: pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                                pb.ETA()])(x)


def get_tcc_values(data, force_new=False):
    if 'tcc' in data.root and not force_new:
        return data.root.tcc.read()
    else:
        h_events = data.root.hisparc.cluster_kascade.station_601.events
        c_index = data.root.kascade.c_index

        tcc = []
        for idx in progressbar(c_index):
            event = h_events[idx['h_idx']]
            value = calculate_tcc(event)
            tcc.append(value)

        tcc = array(tcc)
        if 'tcc' in data.root:
            data.removeNode('/tcc')
        data.createArray('/', 'tcc', tcc)
        return tcc


def calculate_tcc(event):
    n = array([event[u] for u in 'n1', 'n2', 'n3', 'n4'])
    n = where(n < .5, 0, n)
    if not (n >= 0).all():
        return -999
    
    i_max = len(n)
    mean_n = n.mean()

    if mean_n == 0.:
        return -998

    variance = ((n - mean_n) ** 2).sum() / (i_max - 1)
    T_S = variance / mean_n
    T_CC = (i_max - 1) * T_S
    return T_CC


def plot_core_positions(data):
    core_dist = _get_core_dists(data, 'tcc >= 10')
    N = len(core_dist)
    false_core_dist = _get_core_dists(data, 'tcc >= 0', limit=N)
    small_tcc_core_dist = _get_core_dists(data, 'tcc < 10', limit=N)

    figure()
    hist(core_dist, bins=50, histtype='step', label='tcc >= 10')
    hist(false_core_dist, bins=50, histtype='step', label='uncorrelated')
    hist(small_tcc_core_dist, bins=50, histtype='step', label='tcc < 10')
    legend()
    xlabel("Core distance [m]")
    ylabel("Counts")

    graph = artist.GraphArtist()
    n, bins = histogram(core_dist, bins=linspace(0, 200, 51))
    graph.histogram(n, bins)
    graph.add_pin(r'$T_{CC} \geq 10$', x=18, location='above right',
                  use_arrow=True)
    n, bins = histogram(false_core_dist, bins=linspace(0, 200, 51))
    graph.histogram(n, bins, linestyle='gray')
    graph.add_pin('uncorrelated', x=37, location='above right',
                  use_arrow=True)
    graph.set_xlabel(r"Core distance [\si{\meter}]")
    graph.set_ylabel("Counts")
    graph.set_xlimits(0, 200)
    graph.set_ylimits(min=0)
    graph.save_as_pdf('preview')


def scatter_core_positions(data):
    cores = _get_core_positions(data, 'tcc >= 10')
    N = len(cores)
    small_tcc_cores = _get_core_positions(data, 'tcc < 10', limit=N)

    figure()
    x, y = zip(*cores)
    plot(x, y, ',')

    figure()
    x, y = zip(*small_tcc_cores)
    plot(x, y, ',')


def _get_core_dists(data, sel_str, limit=None):

    core_pos = _get_core_positions(data, sel_str, limit)
    x, y = zip(*core_pos)
    core_dist = sqrt((array(x) - X) ** 2 + (array(y) - Y) ** 2)

    return core_dist


def _get_core_positions(data, sel_str, limit=None):
    h_events = data.root.hisparc.cluster_kascade.station_601.events
    k_events = data.root.kascade.events
    c_index = data.root.kascade.c_index
    tcc = data.root.tcc.read()

    sel = eval(sel_str)
    sel_indices = sel.nonzero()[0]
    indices = c_index.readCoordinates(sel_indices)
    k_idx = indices['k_idx']

    core_pos = k_events.readCoordinates(k_idx, field='core_pos')
    if limit:
        core_pos = core_pos[:limit]

    return core_pos


def plot_energy(data, sel_str):
    h_events = data.root.hisparc.cluster_kascade.station_601.events
    k_events = data.root.kascade.events
    c_index = data.root.kascade.c_index
    tcc = data.root.tcc.read()

    sel = eval(sel_str)
    sel_indices = sel.nonzero()[0]
    indices = c_index.readCoordinates(sel_indices)
    k_idx = indices['k_idx']

    energy = k_events.readCoordinates(k_idx, field='energy')
    false_energy = k_events.readCoordinates(k_idx + 1, field='energy')

    print len(energy), len(false_energy)

    figure()
    hist(log10(energy), bins=linspace(14, 18, 51), histtype='step', label=sel_str)
    hist(log10(false_energy), bins=linspace(14, 18, 51), histtype='step', label='uncorrelated')
    legend()

    graph = artist.GraphArtist()
    n, bins = histogram(log10(energy), bins=linspace(14, 18, 51))
    graph.histogram(n, bins)
    graph.add_pin(r'$T_{CC} \geq 10$', x=14.6, location='above right',
                  use_arrow=True)
    n, bins = histogram(log10(false_energy), bins=linspace(14, 18, 51))
    graph.histogram(n, bins, linestyle='gray')
    graph.add_pin('uncorrelated', x=15.5, location='above right',
                  use_arrow=True)
    graph.set_xlabel(r"$\lg$ energy [$\lg\si{\electronvolt}$]")
    graph.set_ylabel("Counts")
    graph.set_xlimits(14, 18)
    graph.set_ylimits(min=0)
    graph.save_as_pdf('preview')


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.openFile('kascade.h5', 'a')

    tcc = get_tcc_values(data, force_new=False)
    #plot_core_positions(data)
    #scatter_core_positions(data)
    plot_energy(data, 'tcc >= 10')
    #plot_energy(data, 'tcc < 10')
