import itertools

import tables
from pylab import *

import master
import utils


USE_TEX = False

# For matplotlib plots
if USE_TEX:
    rcParams['font.serif'] = 'Computer Modern'
    rcParams['font.sans-serif'] = 'Computer Modern'
    rcParams['font.family'] = 'sans-serif'
    rcParams['figure.figsize'] = [4 * x for x in (1, 2. / 3)]
    rcParams['figure.subplot.left'] = 0.175
    rcParams['figure.subplot.bottom'] = 0.175
    rcParams['font.size'] = 10
    rcParams['legend.fontsize'] = 'small'
    rcParams['text.usetex'] = True


def main(data):
    plot_sciencepark_cluster()
    plot_all_single_and_cluster_combinations(data)

def plot_sciencepark_cluster():
    cluster = master.ScienceParkCluster(range(501, 507))

    figure()
    xl, yl = [], []
    for station in cluster.stations:
        for detector in station.detectors:
            x, y = detector.get_xy_coordinates()
            xl.append(x)
            yl.append(y)
            scatter(x, y, c='black', s=3)
    axis('equal')

    utils.savedata([xl, yl])

def plot_all_single_and_cluster_combinations(data):
    for station_group in itertools.combinations(range(501, 507), 3):
        for station in station_group:
            plot_direction_single_vs_cluster(data, station, station_group)

def plot_direction_single_vs_cluster(data, station, cluster):
    reconstructions = data.root.reconstructions.reconstructions

    station_query = '(N == 1) & s%d' % station
    cluster_query = '(N == 3) & ' + ' & '.join(['s%d' % u for u in cluster])
    cluster_str = [str(u) for u in cluster]

    sel_station = reconstructions.readWhere(station_query)
    sel_cluster = reconstructions.readWhere(cluster_query)

    theta_cluster, theta_station = [], []
    phi_cluster, phi_station = [], []
    for event_cluster in sel_cluster:
        coinc_id = event_cluster['coinc_id']
        event_station = sel_station.compress(sel_station[:]['coinc_id'] == coinc_id)
        assert len(event_station) <= 1
        if event_station:
            theta_cluster.append(event_cluster['reconstructed_theta'])
            phi_cluster.append(event_cluster['reconstructed_phi'])
            theta_station.append(event_station['reconstructed_theta'])
            phi_station.append(event_station['reconstructed_phi'])

    figsize = list(rcParams['figure.figsize'])
    figsize[1] = figsize[0] / 2
    figure(figsize=figsize)
    subplot(121)
    plot(theta_station, theta_cluster, ',')
    xlabel(r"$\theta_{%d}$" % station)
    ylabel(r"$\theta_{\{%s\}}$" % ','.join(cluster_str))
    xlim(0, pi / 2)
    ylim(0, pi / 2)

    subplot(122)
    plot(phi_station, phi_cluster, ',')
    xlabel(r"$\phi_{%d}$" % station)
    ylabel(r"$\phi_{\{%s\}}$" % ','.join(cluster_str))
    xlim(-pi, pi)
    ylim(-pi, pi)

    utils.saveplot('%d-%s' % (station, '_'.join(cluster_str)))


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.openFile('master.h5')

    utils.set_prefix("SP-DIR-")
    main(data)
