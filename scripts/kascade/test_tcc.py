import tables
import time

from sapphire.analysis.core_reconstruction import CoreReconstruction, \
                                                  PlotCoreReconstruction
from sapphire.utils import pbar

import artist


X, Y = 65., 20.82


def get_tcc_values(data, force_new=False):
    if 'tcc' in data.root and not force_new:
        return data.root.tcc.read()
    else:
        h_events = data.root.hisparc.cluster_kascade.station_601.events
        c_index = data.root.kascade.c_index

        tcc = []
        for idx in pbar(c_index):
            event = h_events[idx['h_idx']]
            value = calculate_tcc(event)
            tcc.append(value)

        tcc = array(tcc)
        if 'tcc' in data.root:
            data.remove_node('/tcc')
        data.create_array('/', 'tcc', tcc)
        return tcc


def calculate_tcc(event):
    n = array([event[u] for u in 'n1', 'n2', 'n3', 'n4'])
    n = where(n < .5, 0, n)
    if not (n > 0).sum() >= 2:
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
    indices = c_index.read_coordinates(sel_indices)
    k_idx = indices['k_idx']

    core_pos = k_events.read_coordinates(k_idx, field='core_pos')
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
    indices = c_index.read_coordinates(sel_indices)
    k_idx = indices['k_idx']

    energy = k_events.read_coordinates(k_idx, field='energy')
    false_energy = k_events.read_coordinates(k_idx + 1, field='energy')

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


def reconstruct_shower_sizes(data, tcc):
    reconstruction = KascadeCoreReconstruction(data, '/core',
                                               overwrite=True)
    reconstruction.reconstruct_core_positions(
        '/hisparc/cluster_kascade/station_601', '/kascade', tcc)


class KascadeCoreReconstruction(CoreReconstruction):
    def reconstruct_core_positions(self, hisparc_group, kascade_group, tcc):
        hisparc_group = self.data.get_node(hisparc_group)

        hisparc_table = self.data.get_node(hisparc_group, 'events')
        c_index = self.data.get_node(kascade_group, 'c_index')
        kascade_table = self.data.get_node(kascade_group, 'events')

        self.cluster = hisparc_group._v_attrs.cluster
        self._store_cluster_with_results()

        for idx, tcc_value in pbar(zip(c_index[:self.N], tcc)):
            hisparc_event = hisparc_table[idx['h_idx']]
            kascade_event = kascade_table[idx['k_idx']]

            if tcc_value >= 10:
                x, y, N = self.reconstruct_core_position(hisparc_event)
                self.store_reconstructed_event(hisparc_event,
                                               kascade_event, x, y, N)

        self.results_table.flush()

    def store_reconstructed_event(self, hisparc_event, kascade_event,
                                  reconstructed_core_x,
                                  reconstructed_core_y,
                                  reconstructed_shower_size):
        dst_row = self.results_table.row

        dst_row['id'] = hisparc_event['event_id']
        dst_row['station_id'] = 0
        dst_row['t1'] = hisparc_event['t1']
        dst_row['t2'] = hisparc_event['t2']
        dst_row['t3'] = hisparc_event['t3']
        dst_row['t4'] = hisparc_event['t4']
        dst_row['n1'] = hisparc_event['n1']
        dst_row['n2'] = hisparc_event['n2']
        dst_row['n3'] = hisparc_event['n3']
        dst_row['n4'] = hisparc_event['n4']
        dst_row['reference_theta'] = kascade_event['zenith']
        dst_row['reference_phi'] = kascade_event['azimuth']
        dst_row['reference_core_pos'] = kascade_event['core_pos']
        dst_row['reconstructed_core_pos'] = reconstructed_core_x, \
                                            reconstructed_core_y
        dst_row['reference_shower_size'] = kascade_event['Num_e']
        dst_row['reconstructed_shower_size'] = reconstructed_shower_size
        dst_row['min_n134'] = min(hisparc_event['n1'],
                                  hisparc_event['n3'],
                                  hisparc_event['n4'])
        dst_row.append()

    def get_events_from_coincidence(self, coincidence):
        """Fake coincidences"""
        events = [coincidence]
        return events

    def _get_station_from_event(self, event):
        return self.cluster.stations[0]

    def _station_has_triggered(self, event):
        return True


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.open_file('kascade.h5', 'a')

    tcc = get_tcc_values(data, force_new=False)
    reconstruct_shower_sizes(data, tcc)
    core = data.root.core

    #plot_core_positions(data)
    #scatter_core_positions(data)
    #plot_energy(data, 'tcc >= 10')
    #plot_energy(data, 'tcc < 10')
