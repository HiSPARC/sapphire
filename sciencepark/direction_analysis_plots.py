import itertools

import tables
from pylab import *
from scipy.stats import scoreatpercentile

from sapphire import clusters
from sapphire.analysis.direction_reconstruction import DirectionReconstruction
import utils


STATION_TIMING_ERR = 3.1
CLUSTER_TIMING_ERR = 5.9


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
    #plot_sciencepark_cluster()
    #plot_all_single_and_cluster_combinations(data)
    #hist_phi_single_stations(data)
    #hist_theta_single_stations(data)
    #plot_N_vs_R(data)
    #plot_fav_single_vs_cluster(data)
    plot_fav_uncertainty_single_vs_cluster(data)
    #hist_fav_single_stations(data)

def plot_sciencepark_cluster():
    cluster = clusters.ScienceParkCluster(range(501, 507))

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
    utils.saveplot()

def plot_all_single_and_cluster_combinations(data):
    for station_group in itertools.combinations(range(501, 507), 3):
        for station in station_group:
            plot_direction_single_vs_cluster(data, station, station_group)

def calc_direction_single_vs_cluster(data, station, cluster):
    reconstructions = data.root.reconstructions.reconstructions

    station_query = '(N == 1) & s%d & (min_n134 >= 2.)' % station
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

    return array(theta_station).flatten(), array(phi_station).flatten(), \
        array(theta_cluster).flatten(), array(phi_cluster).flatten()

def plot_direction_single_vs_cluster(data, station, cluster):
    cluster_str = [str(u) for u in cluster]

    theta_station, phi_station, theta_cluster, phi_cluster = \
        calc_direction_single_vs_cluster(data, station, cluster)

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

def hist_phi_single_stations(data):
    reconstructions = data.root.reconstructions.reconstructions

    figure()
    for n, station in enumerate(range(501, 507), 1):
        subplot(2, 3, n)
        query = '(N == 1) & s%d' % station
        phi = reconstructions.readWhere(query, field='reconstructed_phi')
        hist(rad2deg(phi), bins=linspace(-180, 180, 21), histtype='step')
        xlabel(r"$\phi$")
        legend([station])
        locator_params(tight=True, nbins=4)

    utils.saveplot()

def hist_theta_single_stations(data):
    reconstructions = data.root.reconstructions.reconstructions

    figure()
    for n, station in enumerate(range(501, 507), 1):
        subplot(2, 3, n)
        query = '(N == 1) & s%d' % station
        theta = reconstructions.readWhere(query, field='reconstructed_theta')
        hist(rad2deg(theta), bins=linspace(0, 45, 21), histtype='step')
        xlabel(r"$\theta$")
        legend([station])
        locator_params(tight=True, nbins=4)

    utils.saveplot()

def plot_N_vs_R(data):
    stations = range(501, 507)
    station_ids = range(6)
    cluster = clusters.ScienceParkCluster(stations)

    c_index = data.root.coincidences.c_index
    observables = data.root.coincidences.observables

    stations_in_coincidence = []
    for coincidence_events in c_index:
        stations = [observables[u]['station_id'] for u in
                    coincidence_events]
        stations_in_coincidence.append(stations)

    figure()
    for station1, station2 in itertools.combinations(station_ids, 2):
        condition = [station1 in u and station2 in u for u in
                     stations_in_coincidence]
        N = sum(condition)
        R, phi = cluster.calc_r_and_phi_for_stations(station1, station2)
        scatter(R, N)
    xlabel("Distance [m]")
    ylabel("Number of coincidences")

    utils.saveplot()

def plot_fav_single_vs_cluster(data):
    cluster = [501, 503, 506]
    cluster_str = [str(u) for u in cluster]

    figure()
    for n, station in enumerate(cluster, 1):
        theta_station, phi_station, theta_cluster, phi_cluster = \
            calc_direction_single_vs_cluster(data, station, cluster)

        subplot(2, 3, n)
        plot(rad2deg(phi_station), rad2deg(phi_cluster), ',')
        xlabel(r"$\phi_{%d}$" % station)
        ylabel(r"$\phi_{\{%s\}}$" % ','.join(cluster_str))
        xlim(-180, 180)
        ylim(-180, 180)
        locator_params(tight=True, nbins=4)

        subplot(2, 3, n + 3)
        plot(rad2deg(theta_station), rad2deg(theta_cluster), ',')
        xlabel(r"$\theta_{%d}$" % station)
        ylabel(r"$\theta_{\{%s\}}$" % ','.join(cluster_str))
        xlim(0, 45)
        ylim(0, 45)
        locator_params(tight=True, nbins=4)
    utils.saveplot()

def plot_fav_uncertainty_single_vs_cluster(data):
    cluster = [501, 503, 506]
    cluster_str = [str(u) for u in cluster]
    cluster_ids = [0, 2, 5]

    figure()
    for n, station in enumerate(cluster, 1):
        theta_station, phi_station, theta_cluster, phi_cluster = \
            calc_direction_single_vs_cluster(data, station, cluster)

        bins = linspace(0, deg2rad(35), 11)
        x, y, y2 = [], [], []
        for low, high in zip(bins[:-1], bins[1:]):
            sel_phi_c = phi_cluster.compress((low <= theta_station) &
                                             (theta_station < high))
            sel_phi_s = phi_station.compress((low <= theta_station) &
                                             (theta_station < high))
            sel_theta_c = theta_cluster.compress((low <= theta_station) &
                                                 (theta_station < high))
            sel_theta_s = theta_station.compress((low <= theta_station) &
                                                 (theta_station < high))
            dphi = sel_phi_s - sel_phi_c
            dtheta = sel_theta_s - sel_theta_c
            # make sure phi, theta are between -pi and pi
            dphi = (dphi + pi) % (2 * pi) - pi
            dtheta = (dtheta + pi) % (2 * pi) - pi
            print rad2deg((low + high) / 2), len(dphi), len(dtheta)
            x.append((low + high) / 2)
            #y.append(std(dphi))
            #y2.append(std(dtheta))
            y.append((scoreatpercentile(dphi, 83) - scoreatpercentile(dphi, 17)) / 2)
            y2.append((scoreatpercentile(dtheta, 83) - scoreatpercentile(dtheta, 17)) / 2)

        ex = linspace(0, deg2rad(35), 50)
        ephi, etheta = [], []
        for theta in ex:
            ephi.append(calc_phi_error_for_station_cluster(theta, n,
                                                           cluster_ids))
            etheta.append(calc_theta_error_for_station_cluster(theta, n,
                                                               cluster_ids))

        subplot(2, 3, n)
        plot(rad2deg(x), rad2deg(y))
        plot(rad2deg(ex), rad2deg(ephi))
        xlabel(r"$\theta_{%d}$ [deg]" % station)
        if n == 1:
            ylabel(r"$\phi$ uncertainty [deg]")
        ylim(0, 100)
        locator_params(tight=True, nbins=4)

        subplot(2, 3, n + 3)
        plot(rad2deg(x), rad2deg(y2))
        plot(rad2deg(ex), rad2deg(etheta))
        xlabel(r"$\theta_{%d}$ [deg]" % station)
        if n == 1:
            ylabel(r"$\theta$ uncertainty [deg]")
        #ylabel(r"$\theta_{\{%s\}}$" % ','.join(cluster_str))
        ylim(0, 15)
        locator_params(tight=True, nbins=4)
    utils.saveplot()

def calc_phi_error_for_station_cluster(theta, station, cluster):
    phis = linspace(-pi, pi, 50)
    rec = DirectionReconstruction
    sciencepark = clusters.ScienceParkCluster(range(501, 507))

    r1, phi1 = sciencepark.stations[station].calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = sciencepark.stations[station].calc_r_and_phi_for_detectors(1, 4)
    err_single = rec.rel_phi_errorsq(theta, phis, phi1, phi2, r1, r2)

    r1, phi1 = sciencepark.calc_r_and_phi_for_stations(cluster[0], cluster[1])
    r2, phi2 = sciencepark.calc_r_and_phi_for_stations(cluster[0], cluster[2])
    err_cluster = rec.rel_phi_errorsq(theta, phis, phi1, phi2, r1, r2)

    # errors are already squared!!
    err_total = sqrt(STATION_TIMING_ERR ** 2 * err_single +
                     CLUSTER_TIMING_ERR ** 2 * err_cluster)
    return mean(err_total)

def calc_theta_error_for_station_cluster(theta, station, cluster):
    phis = linspace(-pi, pi, 50)
    rec = DirectionReconstruction
    sciencepark = clusters.ScienceParkCluster(range(501, 507))

    r1, phi1 = sciencepark.stations[station].calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = sciencepark.stations[station].calc_r_and_phi_for_detectors(1, 4)
    err_single = rec.rel_theta1_errorsq(theta, phis, phi1, phi2, r1, r2)

    r1, phi1 = sciencepark.calc_r_and_phi_for_stations(cluster[0], cluster[1])
    r2, phi2 = sciencepark.calc_r_and_phi_for_stations(cluster[0], cluster[2])
    err_cluster = rec.rel_theta1_errorsq(theta, phis, phi1, phi2, r1, r2)

    # errors are already squared!!
    err_total = sqrt(STATION_TIMING_ERR ** 2 * err_single +
                     CLUSTER_TIMING_ERR ** 2 * err_cluster)
    return mean(err_total)

def hist_fav_single_stations(data):
    reconstructions = data.root.reconstructions.reconstructions

    figure()
    for n, station in enumerate([501, 503, 506], 1):
        query = '(N == 1) & s%d' % station
        phi = reconstructions.readWhere(query, field='reconstructed_phi')
        theta = reconstructions.readWhere(query, field='reconstructed_theta')

        subplot(2, 3, n)
        hist(rad2deg(phi), bins=linspace(-180, 180, 21), histtype='step')
        xlabel(r"$\phi$")
        legend([station], loc='lower right')
        locator_params(tight=True, nbins=4)

        subplot(2, 3, n + 3)
        hist(rad2deg(theta), bins=linspace(0, 45, 21), histtype='step')
        xlabel(r"$\theta$")
        legend([station], loc='lower right')
        locator_params(tight=True, nbins=4)

    utils.saveplot()


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.openFile('new.h5')

    utils.set_prefix("SP-DIR-")
    main(data)
