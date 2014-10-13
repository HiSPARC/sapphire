from __future__ import division

import tables
from itertools import izip
import os.path

from pylab import *

from scipy import integrate
from scipy.special import erf
from scipy.stats import scoreatpercentile
from scipy.interpolate import spline

import utils

from sapphire.analysis.direction_reconstruction import (DirectionReconstruction,
                                                        BinnedDirectionReconstruction)
from sapphire.utils import pbar
from myshowerfront import *

from artist import GraphArtist
import artist.utils


USE_TEX = False

ONEP_TIMING_ERROR = 4.6
TIMING_ERROR = 1.8


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


def do_full_reconstruction(data, N=None):
    """Do a reconstruction of all simulated data and store results"""

    if '/reconstructions' in data:
        data.remove_node('/reconstructions', recursive=True)

    simulation_group = '/simulations/E_1PeV'
    reconstruction_group = simulation_group.replace('simulations', 'reconstructions')

    for theta in 0, 5, 10, 15, 22.5, 30, 35, 45:
        tablename = 'zenith_%s' % str(theta).replace('.', '_')
        source = os.path.join(simulation_group, tablename)
        dest = os.path.join(reconstruction_group, tablename)

        rec = DirectionReconstruction(data, dest, min_n134=1, N=N)
        rec.reconstruct_angles_for_shower_group(source)

    # SPECIALS
    # Station sizes
    for size in 5, 20:
        tablename = 'zenith_22_5_size%d' % size
        source = os.path.join(simulation_group, tablename)
        dest = os.path.join(reconstruction_group, tablename)

        rec = DirectionReconstruction(data, dest, min_n134=1, N=N)
        rec.reconstruct_angles_for_shower_group(source)

    # SPECIALS
    # Binnings

    source = os.path.join(simulation_group, 'zenith_22_5')
    for randomize in False, True:
        dest = os.path.join(reconstruction_group, 'zenith_22_5_binned')
        if randomize:
            dest += '_randomized'

        for binning in 1, 2.5, 5:
            dest_table = dest + '_%s' % str(binning).replace('.', '_')
            rec = BinnedDirectionReconstruction(data, dest_table, min_n134=1, binning=binning, randomize_binning=randomize, N=N)
            rec.reconstruct_angles_for_shower_group(source)


def do_reconstruction_plots(data):
    """Make plots based upon earlier reconstructions"""

    group = data.root.reconstructions

    plot_uncertainty_mip(group)
    plot_uncertainty_zenith(group)
    plot_uncertainty_core_distance(group)
    plot_uncertainty_size(group)
    plot_uncertainty_binsize(group)
    plot_uncertainty_zenith_angular_distance(group)

    plot_phi_reconstruction_results_for_MIP(group, 1)
    plot_phi_reconstruction_results_for_MIP(group, 2)
    boxplot_theta_reconstruction_results_for_MIP(group, 1)
    boxplot_theta_reconstruction_results_for_MIP(group, 2)
    boxplot_phi_reconstruction_results_for_MIP(group, 1)
    boxplot_phi_reconstruction_results_for_MIP(group, 2)
    boxplot_arrival_times(group, 1)
    boxplot_arrival_times(group, 2)
    boxplot_core_distances_for_mips(group)
    save_for_kascade_boxplot_core_distances_for_mips(group)
    plot_detection_efficiency_vs_R_for_angles(1)
    plot_detection_efficiency_vs_R_for_angles(2)
    #plot_reconstruction_efficiency_vs_R_for_angles(1)
    #plot_reconstruction_efficiency_vs_R_for_angles(2)
    artistplot_reconstruction_efficiency_vs_R_for_angles(1)
    artistplot_reconstruction_efficiency_vs_R_for_angles(2)
    #plot_reconstruction_efficiency_vs_R_for_mips()


def plot_uncertainty_mip(group):
    table = group.E_1PeV.zenith_22_5
    rec = DirectionReconstruction

    # constants for uncertainty estimation
    station = table.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    R_list = get_median_core_distances_for_mips(group, range(1, 6))

    figure()
    x, y, y2 = [], [], []
    for N in range(1, 5):
        x.append(N)
        events = table.read_where('min_n134>=%d' % N)
        #query = '(n1 == N) & (n3 == N) & (n4 == N)'
        #vents = table.read_where(query)
        print len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        #y.append(std(errors))
        #y2.append(std(errors2))
        y.append((scoreatpercentile(errors, 83) - scoreatpercentile(errors, 17)) / 2)
        y2.append((scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17)) / 2)
        print "YYY", rad2deg(scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17))

    plot(x, rad2deg(y), '^', label="Theta")
    plot(x, rad2deg(y2), 'v', label="Phi")
    Sx = x
    Sy = y
    Sy2 = y2
    print
    print "mip: min_n134, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    utils.savedata((x, y, y2))

    # Uncertainty estimate
    x = [1, 2, 3, 4, 5]
    phis = linspace(-pi, pi, 50)
    phi_errsq = mean(rec.rel_phi_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    theta_errsq = mean(rec.rel_theta1_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    y = TIMING_ERROR * std_t(x) * sqrt(phi_errsq)
    y2 = TIMING_ERROR * std_t(x) * sqrt(theta_errsq)

    mc = my_std_t_for_R(data, x, R_list)
    for u, v in zip(mc, R_list):
        print v, u, sqrt(u ** 2 + 1.2 ** 2), sqrt((.66 * u) ** 2 + 1.2 ** 2)
    mc = sqrt(mc ** 2 + 1.2 ** 2)
    y3 = mc * sqrt(phi_errsq)
    y4 = mc * sqrt(theta_errsq)

    nx = linspace(1, 4, 100)
    y = spline(x, y, nx)
    y2 = spline(x, y2, nx)
    y3 = spline(x, y3, nx)
    y4 = spline(x, y4, nx)

    plot(nx, rad2deg(y), label="Gauss Phi")
    plot(nx, rad2deg(y2), label="Gauss Theta")
    plot(nx, rad2deg(y3), label="Monte Carlo Phi")
    plot(nx, rad2deg(y4), label="Monte Carlo Theta")
    # Labels etc.
    xlabel("Minimum number of particles")
    ylabel("Angle reconstruction uncertainty [deg]")
    #title(r"$\theta = 22.5^\circ$")
    legend(numpoints=1)
    xlim(.5, 4.5)
    utils.saveplot()
    print

    graph = GraphArtist()
    graph.plot(Sx, rad2deg(Sy), mark='o', linestyle='only marks')
    graph.plot(Sx, rad2deg(Sy2), mark='*', linestyle='only marks')
    graph.plot(nx, rad2deg(y), mark=None, linestyle='dashed,smooth')
    graph.plot(nx, rad2deg(y2), mark=None, linestyle='dashed,smooth')
    graph.set_xlabel("Minimum number of particles")
    graph.set_ylabel(r"Reconstruction uncertainty [\si{\degree}]")
    graph.set_xticks(range(1, 5))
    graph.set_ylimits(0, 32)
    artist.utils.save_graph(graph, dirname='plots')
    graph.plot(nx, rad2deg(y3), mark=None, linestyle='smooth')
    graph.plot(nx, rad2deg(y4), mark=None, linestyle='smooth')
    artist.utils.save_graph(graph, suffix='full', dirname='plots')


def plot_uncertainty_zenith(group):
    group = group.E_1PeV
    rec = DirectionReconstruction

    N = 2

    graph = GraphArtist()

    # constants for uncertainty estimation
    # BEWARE: stations must be the same over all reconstruction tables used
    station = group.zenith_0.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    figure()
    x, y, y2 = [], [], []
    for THETA in 0, 5, 10, 15, 22.5, 30, 35, 45:
        x.append(THETA)
        table = group._f_get_child('zenith_%s' % str(THETA).replace('.', '_'))
        events = table.read_where('min_n134 >= N')
        print THETA, len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        #y.append(std(errors))
        #y2.append(std(errors2))
        y.append((scoreatpercentile(errors, 83) - scoreatpercentile(errors, 17)) / 2)
        y2.append((scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17)) / 2)
    plot(x, rad2deg(y), '^', label="Theta")
    graph.plot(x, rad2deg(y), mark='o', linestyle=None)
    # Azimuthal angle undefined for zenith = 0
    plot(x[1:], rad2deg(y2[1:]), 'v', label="Phi")
    graph.plot(x[1:], rad2deg(y2[1:]), mark='*', linestyle=None)
    print
    print "zenith: theta, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    utils.savedata((x, y, y2))

    # Uncertainty estimate
    x = linspace(0, deg2rad(45), 50)
    phis = linspace(-pi, pi, 50)
    y, y2 = [], []
    for t in x:
        y.append(mean(rec.rel_phi_errorsq(t, phis, phi1, phi2, r1, r2)))
        y2.append(mean(rec.rel_theta1_errorsq(t, phis, phi1, phi2, r1, r2)))
    y = TIMING_ERROR * sqrt(array(y))
    y2 = TIMING_ERROR * sqrt(array(y2))
    plot(rad2deg(x), rad2deg(y), label="Estimate Phi")
    graph.plot(rad2deg(x), rad2deg(y), mark=None)
    plot(rad2deg(x), rad2deg(y2), label="Estimate Theta")
    graph.plot(rad2deg(x), rad2deg(y2), mark=None)

    # Labels etc.
    xlabel("Shower zenith angle [deg]")
    graph.set_xlabel(r"Shower zenith angle [\si{\degree}]")
    ylabel("Angle reconstruction uncertainty [deg]")
    graph.set_ylabel(r"Angle reconstruction uncertainty [\si{\degree}]")
    #title(r"$N_{MIP} \geq %d$" % N)
    ylim(0, 100)
    graph.set_ylimits(0, 60)
    legend(numpoints=1)
    utils.saveplot()
    artist.utils.save_graph(graph, dirname='plots')
    print


def plot_uncertainty_core_distance(group):
    table = group.E_1PeV.zenith_22_5

    N = 2
    DR = 10

    figure()
    x, y, y2 = [], [], []
    for R in range(0, 81, 20):
        x.append(R)
        events = table.read_where('(min_n134 == N) & (abs(r - R) <= DR)')
        print len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        #y.append(std(errors))
        #y2.append(std(errors2))
        y.append((scoreatpercentile(errors, 83) - scoreatpercentile(errors, 17)) / 2)
        y2.append((scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17)) / 2)

    print
    print "R: theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    utils.savedata((x, y, y2))

    # Plots
    plot(x, rad2deg(y), '^-', label="Theta")
    plot(x, rad2deg(y2), 'v-', label="Phi")

    # Labels etc.
    xlabel("Core distance [m] $\pm %d$" % DR)
    ylabel("Angle reconstruction uncertainty [deg]")
    #title(r"$N_{MIP} = %d, \theta = 22.5^\circ$" % N)
    ylim(ymin=0)
    legend(numpoints=1, loc='best')
    utils.saveplot()
    print


def plot_uncertainty_size(group):
    group = group.E_1PeV
    rec = DirectionReconstruction

    N = 2

    graph = GraphArtist()

    # constants for uncertainty estimation
    # BEWARE: stations must be the same shape(!) over all reconstruction tables used
    station = group.zenith_22_5.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)
    del r1, r2

    figure()
    x, y, y2 = [], [], []
    for size in [5, 10, 20]:
        x.append(size)
        if size != 10:
            table = group._f_get_child('zenith_22_5_size%d' % size)
        else:
            table = group._f_get_child('zenith_22_5')

        events = table.read_where('min_n134 >= N')
        print size, len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        #y.append(std(errors))
        #y2.append(std(errors2))
        y.append((scoreatpercentile(errors, 83) - scoreatpercentile(errors, 17)) / 2)
        y2.append((scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17)) / 2)
    plot(x, rad2deg(y), '^', label="Theta")
    graph.plot(x, rad2deg(y), mark='o', linestyle=None)
    plot(x, rad2deg(y2), 'v', label="Phi")
    graph.plot(x, rad2deg(y2), mark='*', linestyle=None)
    print
    print "stationsize: size, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print

    # Uncertainty estimate
    x = linspace(5, 20, 50)
    phis = linspace(-pi, pi, 50)
    y, y2 = [], []
    for s in x:
        y.append(mean(rec.rel_phi_errorsq(pi / 8, phis, phi1, phi2, r1=s, r2=s)))
        y2.append(mean(rec.rel_theta1_errorsq(pi / 8, phis, phi1, phi2, r1=s, r2=s)))
    y = TIMING_ERROR * sqrt(array(y))
    y2 = TIMING_ERROR * sqrt(array(y2))
    plot(x, rad2deg(y), label="Estimate Phi")
    graph.plot(x, rad2deg(y), mark=None)
    plot(x, rad2deg(y2), label="Estimate Theta")
    graph.plot(x, rad2deg(y2), mark=None)

    # Labels etc.
    xlabel("Station size [m]")
    graph.set_xlabel(r"Station size [\si{\meter}]")
    ylabel("Angle reconstruction uncertainty [deg]")
    graph.set_ylabel(r"Angle reconstruction uncertainty [\si{\degree}]")
    graph.set_ylimits(0, 25)
    #title(r"$\theta = 22.5^\circ, N_{MIP} \geq %d$" % N)
    legend(numpoints=1)
    utils.saveplot()
    artist.utils.save_graph(graph, dirname='plots')
    print


def plot_uncertainty_binsize(group):
    group = group.E_1PeV
    rec = DirectionReconstruction

    N = 2

    graph = GraphArtist()

    # constants for uncertainty estimation
    # BEWARE: stations must be the same over all reconstruction tables used
    station = group.zenith_22_5.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    figure()
    x, y, y2 = [], [], []
    for bin_size in [0, 1, 2.5, 5]:
        x.append(bin_size)
        if bin_size != 0:
            table = group._f_get_child('zenith_22_5_binned_randomized_%s' % str(bin_size).replace('.', '_'))
        else:
            table = group.zenith_22_5
        events = table.read_where('min_n134 >= 2')

        print bin_size, len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        #y.append(std(errors))
        #y2.append(std(errors2))
        y.append((scoreatpercentile(errors, 83) - scoreatpercentile(errors, 17)) / 2)
        y2.append((scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17)) / 2)
    plot(x, rad2deg(y), '^', label="Theta")
    graph.plot(x, rad2deg(y), mark='o', linestyle=None)
    plot(x, rad2deg(y2), 'v', label="Phi")
    graph.plot(x, rad2deg(y2), mark='*', linestyle=None)
    print
    print "binsize: size, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print

    # Uncertainty estimate
    x = linspace(0, 5, 50)
    phis = linspace(-pi, pi, 50)
    y, y2 = [], []
    phi_errorsq = mean(rec.rel_phi_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    theta_errorsq = mean(rec.rel_theta1_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    for t in x:
        y.append(sqrt((TIMING_ERROR ** 2 + t ** 2 / 12) * phi_errorsq))
        y2.append(sqrt((TIMING_ERROR ** 2 + t ** 2 / 12) * theta_errorsq))
    y = array(y)
    y2 = array(y2)
    plot(x, rad2deg(y), label="Estimate Phi")
    graph.plot(x, rad2deg(y), mark=None)
    plot(x, rad2deg(y2), label="Estimate Theta")
    graph.plot(x, rad2deg(y2), mark=None)

    # Labels etc.
    xlabel("Sampling time [ns]")
    graph.set_xlabel(r"Sampling time [\si{\nano\second}]")
    ylabel("Angle reconstruction uncertainty [deg]")
    graph.set_ylabel(r"Angle reconstruction uncertainty [\si{\degree}]")
    graph.set_ylimits(0, 20)
    #title(r"$\theta = 22.5^\circ, N_{MIP} \geq %d$" % N)
    legend(loc='upper left', numpoints=1)
    ylim(0, 20)
    xlim(-0.1, 5.5)
    utils.saveplot()
    artist.utils.save_graph(graph, dirname='plots')
    print

# Time of first hit pamflet functions
Q = lambda t, n: ((.5 * (1 - erf(t / sqrt(2)))) ** (n - 1)
                  * exp(-.5 * t ** 2) / sqrt(2 * pi))

expv_t = vectorize(lambda n: integrate.quad(lambda t: t * Q(t, n)
                                                      / n ** -1,
                                            - inf, +inf))
expv_tv = lambda n: expv_t(n)[0]
expv_tsq = vectorize(lambda n: integrate.quad(lambda t: t ** 2 * Q(t, n)
                                                        / n ** -1,
                                              - inf, +inf))
expv_tsqv = lambda n: expv_tsq(n)[0]

std_t = lambda n: sqrt(expv_tsqv(n) - expv_tv(n) ** 2)


def plot_phi_reconstruction_results_for_MIP(group, N):
    table = group.E_1PeV.zenith_22_5

    events = table.read_where('min_n134 >= %d' % N)
    sim_phi = events['reference_phi']
    r_phi = events['reconstructed_phi']

    figure()
    plot_2d_histogram(rad2deg(sim_phi), rad2deg(r_phi), 180)
    xlabel(r"$\phi_{simulated}$ [deg]")
    ylabel(r"$\phi_{reconstructed}$ [deg]")
    #title(r"$N_{MIP} \geq %d, \quad \theta = 22.5^\circ$" % N)
    utils.saveplot(N)

    graph = artist.GraphArtist()
    bins = linspace(-180, 180, 73)
    H, x_edges, y_edges = histogram2d(rad2deg(sim_phi), rad2deg(r_phi),
                                      bins=bins)
    graph.histogram2d(H, x_edges, y_edges, type='reverse_bw')
    graph.set_xlabel(r'$\phi_\mathrm{sim}$ [\si{\degree}]')
    graph.set_ylabel(r'$\phi_\mathrm{rec}$ [\si{\degree}]')
    graph.set_xticks(range(-180, 181, 90))
    graph.set_yticks(range(-180, 181, 90))
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def boxplot_theta_reconstruction_results_for_MIP(group, N):
    group = group.E_1PeV

    figure()

    angles = [0, 5, 10, 15, 22.5, 30, 35, 45]
    r_dtheta = []
    d25, d50, d75 = [], [], []
    for angle in angles:
        table = group._f_get_child('zenith_%s' % str(angle).replace('.', '_'))
        sel = table.read_where('min_n134 >= %d' % N)
        dtheta = sel[:]['reconstructed_theta'] - sel[:]['reference_theta']
        r_dtheta.append(rad2deg(dtheta))

        d25.append(scoreatpercentile(rad2deg(dtheta), 25))
        d50.append(scoreatpercentile(rad2deg(dtheta), 50))
        d75.append(scoreatpercentile(rad2deg(dtheta), 75))

    fill_between(angles, d25, d75, color='0.75')
    plot(angles, d50, 'o-', color='black')

    xlabel(r"$\theta_{simulated}$ [deg]")
    ylabel(r"$\theta_{reconstructed} - \theta_{simulated}$ [deg]")
    #title(r"$N_{MIP} \geq %d$" % N)

    axhline(0, color='black')
    ylim(-10, 25)

    utils.saveplot(N)

    graph = GraphArtist()
    graph.draw_horizontal_line(0, linestyle='gray')
    graph.shade_region(angles, d25, d75)
    graph.plot(angles, d50, linestyle=None)
    graph.set_xlabel(r"$\theta_\mathrm{sim}$ [\si{\degree}]")
    graph.set_ylabel(r"$\theta_\mathrm{rec} - \theta_\mathrm{sim}$ [\si{\degree}]")
    graph.set_title(r"$N_\mathrm{MIP} \geq %d$" % N)
    graph.set_ylimits(-8, 22)
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def boxplot_phi_reconstruction_results_for_MIP(group, N):
    table = group.E_1PeV.zenith_22_5

    figure()

    bin_edges = linspace(-180, 180, 18)
    x, r_dphi = [], []
    d25, d50, d75 = [], [], []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        rad_low = deg2rad(low)
        rad_high = deg2rad(high)
        query = '(min_n134 >= N) & (rad_low < reference_phi) & (reference_phi < rad_high)'
        sel = table.read_where(query)
        dphi = sel[:]['reconstructed_phi'] - sel[:]['reference_phi']
        dphi = (dphi + pi) % (2 * pi) - pi
        r_dphi.append(rad2deg(dphi))

        d25.append(scoreatpercentile(rad2deg(dphi), 25))
        d50.append(scoreatpercentile(rad2deg(dphi), 50))
        d75.append(scoreatpercentile(rad2deg(dphi), 75))
        x.append((low + high) / 2)

    fill_between(x, d25, d75, color='0.75')
    plot(x, d50, 'o-', color='black')

    xlabel(r"$\phi_{simulated}$ [deg]")
    ylabel(r"$\phi_{reconstructed} - \phi_{simulated}$ [deg]")
    #title(r"$N_{MIP} \geq %d, \quad \theta = 22.5^\circ$" % N)

    xticks(linspace(-180, 180, 9))
    axhline(0, color='black')
    ylim(-15, 15)

    utils.saveplot(N)

    graph = GraphArtist()
    graph.draw_horizontal_line(0, linestyle='gray')
    graph.shade_region(x, d25, d75)
    graph.plot(x, d50, linestyle=None)
    graph.set_xlabel(r"$\phi_\mathrm{sim}$ [\si{\degree}]")
    graph.set_ylabel(r"$\phi_\mathrm{rec} - \phi_\mathrm{sim}$ [\si{\degree}]")
    graph.set_title(r"$N_\mathrm{MIP} \geq %d$" % N)
    graph.set_xticks([-180, -90, '...', 180])
    graph.set_xlimits(-180, 180)
    graph.set_ylimits(-17, 17)
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def boxplot_arrival_times(group, N):
    table = group.E_1PeV.zenith_0

    sel = table.read_where('min_n134 >= N')
    t1 = sel[:]['t1']
    t3 = sel[:]['t3']
    t4 = sel[:]['t4']
    ts = concatenate([t1, t3, t4])
    print "Median arrival time delay over all detected events", median(ts)

    figure()

    bin_edges = linspace(0, 100, 11)
    x, arrival_times = [], []
    t25, t50, t75 = [], [], []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        query = '(min_n134 >= N) & (low <= r) & (r < high)'
        sel = table.read_where(query)
        t1 = sel[:]['t1']
        t2 = sel[:]['t2']
        ct1 = t1.compress((t1 > -999) & (t2 > -999))
        ct2 = t2.compress((t1 > -999) & (t2 > -999))
        ts = abs(ct2 - ct1)

        t25.append(scoreatpercentile(ts, 25))
        t50.append(scoreatpercentile(ts, 50))
        t75.append(scoreatpercentile(ts, 75))
        x.append((low + high) / 2)

    fill_between(x, t25, t75, color='0.75')
    plot(x, t50, 'o-', color='black')

    xlabel("Core distance [m]")
    ylabel("Arrival time delay [ns]")
    #title(r"$N_{MIP} \geq %d, \quad \theta = 0^\circ$" % N)

    xticks(arange(0, 100.5, 10))

    utils.savedata((x, t25, t50, t75), N)
    utils.saveplot(N)

    graph = GraphArtist()
    graph.shade_region(x, t25, t75)
    graph.plot(x, t50, linestyle=None)
    graph.set_xlabel(r"Core distance [\si{\meter}]")
    graph.set_ylabel(r"Arrival time difference $|t_2 - t_1|$ [\si{\nano\second}]")
    graph.set_xlimits(0, 100)
    graph.set_ylimits(min=0)
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def get_median_core_distances_for_mips(group, N_list):
    table = group.E_1PeV.zenith_22_5

    r_list = []
    r25, r50, r75 = [], [], []
    x = []
    for N in N_list:
        sel = table.read_where('min_n134 >= N')
        #query = '(n1 == N) & (n3 == N) & (n4 == N)'
        #sel = table.read_where(query)
        r = sel[:]['r']
        r_list.append(r)
        x.append(N)

        r25.append(scoreatpercentile(r, 25))
        r50.append(scoreatpercentile(r, 50))
        r75.append(scoreatpercentile(r, 75))
        print "MIP, median, mean", N, r50[-1], mean(r), std(r) / mean(r)

    return r50


def boxplot_core_distances_for_mips(group):
    table = group.E_1PeV.zenith_22_5

    figure()

    r_list = []
    r25, r50, r75 = [], [], []
    x = []
    for N in range(1, 5):
        sel = table.read_where('min_n134 >= N')
        r = sel[:]['r']
        r_list.append(r)
        x.append(N)

        r25.append(scoreatpercentile(r, 25))
        r50.append(scoreatpercentile(r, 50))
        r75.append(scoreatpercentile(r, 75))

    fill_between(x, r25, r75, color='0.75')
    plot(x, r50, 'o-', color='black')

    xticks(range(1, 5))
    xlabel("Minimum number of particles")
    ylabel("Core distance [m]")
    #title(r"$\theta = 22.5^\circ$")

    utils.saveplot()

    graph = GraphArtist()
    graph.shade_region(x, r25, r75)
    graph.plot(x, r50, linestyle=None)
    graph.set_xlabel("Minimum number of particles")
    graph.set_ylabel(r"Core distance [\si{\meter}]")
    graph.set_ylimits(min=0)
    graph.set_xticks(range(5))
    artist.utils.save_graph(graph, dirname='plots')


def save_for_kascade_boxplot_core_distances_for_mips(group):
    table = group.E_1PeV.zenith_22_5

    r25_list = []
    r50_list = []
    r75_list = []
    x = []
    for N in range(1, 5):
        sel = table.read_where('(min_n134 == N) & (r <= 80)')
        r = sel[:]['r']
        r25_list.append(scoreatpercentile(r, 25))
        r50_list.append(scoreatpercentile(r, 50))
        r75_list.append(scoreatpercentile(r, 75))
        x.append(N)

    utils.savedata((x, r25_list, r50_list, r75_list))


def plot_detection_efficiency_vs_R_for_angles(N):
    figure()
    graph = GraphArtist()
    locations = iter(['right', 'left', 'below left'])
    positions = iter([.18, .14, .15])

    bin_edges = linspace(0, 100, 20)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2.

    for angle in [0, 22.5, 35]:
        angle_str = str(angle).replace('.', '_')
        shower_group = '/simulations/E_1PeV/zenith_%s' % angle_str

        efficiencies = []
        for low, high in zip(bin_edges[:-1], bin_edges[1:]):
            shower_results = []
            for shower in data.list_nodes(shower_group):
                sel_query = '(low <= r) & (r < high)'
                coinc_sel = shower.coincidences.read_where(sel_query)
                ids = coinc_sel['id']
                obs_sel = shower.observables.read_coordinates(ids)
                assert (obs_sel['id'] == ids).all()

                o = obs_sel
                sel = obs_sel.compress((o['n1'] >= N) & (o['n3'] >= N) &
                                       (o['n4'] >= N))
                shower_results.append(len(sel) / len(obs_sel))
            efficiencies.append(mean(shower_results))

        plot(x, efficiencies, label=r'$\theta = %s^\circ$' % angle)
        graph.plot(x, efficiencies, mark=None)
        graph.add_pin(r'\SI{%s}{\degree}' % angle,
                      location=locations.next(), use_arrow=True,
                      relative_position=positions.next())

    xlabel("Core distance [m]")
    graph.set_xlabel(r"Core distance [\si{\meter}]")
    ylabel("Detection efficiency")
    graph.set_ylabel("Detection efficiency")
    #title(r"$N_{MIP} \geq %d$" % N)
    legend()
    graph.set_xlimits(0, 100)
    graph.set_ylimits(0, 1)

    utils.saveplot(N)
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def plot_reconstruction_efficiency_vs_R_for_angles(N):
    group = data.root.reconstructions.E_1PeV

    figure()

    bin_edges = linspace(0, 100, 10)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2.

    all_data = []

    for angle in [0, 22.5, 35]:
        angle_str = str(angle).replace('.', '_')
        shower_group = '/simulations/E_1PeV/zenith_%s' % angle_str
        reconstructions = group._f_get_child('zenith_%s' % angle_str)

        efficiencies = []
        for low, high in zip(bin_edges[:-1], bin_edges[1:]):
            shower_results = []
            for shower in data.list_nodes(shower_group):
                sel_query = '(low <= r) & (r < high)'
                coinc_sel = shower.coincidences.read_where(sel_query)
                ids = coinc_sel['id']
                obs_sel = shower.observables.read_coordinates(ids)
                assert (obs_sel['id'] == ids).all()

                o = obs_sel
                sel = obs_sel.compress((o['n1'] >= N) & (o['n3'] >= N) &
                                       (o['n4'] >= N))
                shower_results.append(len(sel))
            ssel = reconstructions.read_where('(min_n134 >= N) & (low <= r) & (r < high)')
            efficiencies.append(len(ssel) / sum(shower_results))

        all_data.append(efficiencies)
        plot(x, efficiencies, label=r'$\theta = %s^\circ$' % angle)

    xlabel("Core distance [m]")
    ylabel("Reconstruction efficiency")
    #title(r"$N_{MIP} \geq %d$" % N)
    legend()

    utils.saveplot(N)
    utils.savedata(array([x] + all_data).T, suffix=N)


def artistplot_reconstruction_efficiency_vs_R_for_angles(N):
    filename = 'DIR-plot_reconstruction_efficiency_vs_R_for_angles-%d.txt' % N
    all_data = loadtxt(os.path.join('plots/', filename))

    graph = GraphArtist()
    locations = iter(['above right', 'below left', 'below left'])
    positions = iter([.9, .2, .2])

    x = all_data[:, 0]

    for angle, efficiencies in zip([0, 22.5, 35], all_data[:, 1:].T):
        graph.plot(x, efficiencies, mark=None)
        graph.add_pin(r'\SI{%s}{\degree}' % angle, use_arrow=True,
                      location=locations.next(),
                      relative_position=positions.next())

    graph.set_xlabel("Core distance [\si{\meter}]")
    graph.set_ylabel("Reconstruction efficiency")
    graph.set_xlimits(0, 100)
    graph.set_ylimits(max=1)
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def plot_reconstruction_efficiency_vs_R_for_mips():
    reconstructions = data.root.reconstructions.E_1PeV.zenith_22_5

    figure()

    bin_edges = linspace(0, 100, 10)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2.

    for N in range(1, 5):
        shower_group = '/simulations/E_1PeV/zenith_22_5'

        efficiencies = []
        for low, high in zip(bin_edges[:-1], bin_edges[1:]):
            shower_results = []
            for shower in data.list_nodes(shower_group):
                sel_query = '(low <= r) & (r < high)'
                coinc_sel = shower.coincidences.read_where(sel_query)
                ids = coinc_sel['id']
                obs_sel = shower.observables.read_coordinates(ids)
                assert (obs_sel['id'] == ids).all()

                o = obs_sel
                sel = o.compress(amin(array([o['n1'], o['n3'], o['n4']]), 0) == N)

                shower_results.append(len(sel))
            ssel = reconstructions.read_where('(min_n134 == N) & (low <= r) & (r < high)')
            print sum(shower_results), len(ssel), len(ssel) / sum(shower_results)
            efficiencies.append(len(ssel) / sum(shower_results))

        plot(x, efficiencies, label=r'$N_{MIP} = %d$' % N)

    xlabel("Core distance [m]")
    ylabel("Reconstruction efficiency")
    #title(r"$\theta = 22.5^\circ$")
    legend()

    utils.saveplot()


def plot_2d_histogram(x, y, bins):
    H, xedges, yedges = histogram2d(x, y, bins)
    imshow(H.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           origin='lower left', interpolation='lanczos', aspect='auto',
           cmap=cm.Greys)
    colorbar()


def do_jos_plots(data):
    make_datasets_failed_reconstructions_scatter(data)
    plot_failed_and_successful_scatter_plots()
    plot_failed_histograms()


def make_datasets_failed_reconstructions_scatter(data):
    global dt1, dt2, phis_sim, phis_rec
    global gdt1, gdt2, gphis_sim, gphis_rec

    group = data.root.simulations.E_1PeV.zenith_22_5.shower_0
    observables = group.observables
    coincidences = group.coincidences

    dt1, dt2, phis_sim, phis_rec = [], [], [], []
    gdt1, gdt2, gphis_sim, gphis_rec = [], [], [], []

    for event, coincidence in pbar(izip(observables, coincidences),
                                   len(observables)):
        assert event['id'] == coincidence['id']
        if min(event['n1'], event['n3'], event['n4']) >= 1.:
            theta, phi = reconstruct_angle(event, 10.)
            assert not isnan(phi)

            if isnan(theta):
                dt1.append(event['t1'] - event['t3'])
                dt2.append(event['t1'] - event['t4'])
                phis_sim.append(coincidence['phi'])
                phis_rec.append(phi)
            else:
                gdt1.append(event['t1'] - event['t3'])
                gdt2.append(event['t1'] - event['t4'])
                gphis_sim.append(coincidence['phi'])
                gphis_rec.append(phi)


def plot_failed_and_successful_scatter_plots():
    figure(figsize=(20., 11.5))

    subplot(231)
    plot(gdt1, rad2deg(gphis_sim), ',', c='green')
    plot(dt1, rad2deg(phis_sim), ',', c='red')
    xlabel(r"$t_1 - t_3$ [ns]")
    ylabel(r"$\phi_{sim}$")
    xlim(-200, 200)

    subplot(232)
    plot(gdt2, rad2deg(gphis_sim), ',', c='green')
    plot(dt2, rad2deg(phis_sim), ',', c='red')
    xlabel(r"$t_1 - t_4$ [ns]")
    ylabel(r"$\phi_{sim}$")
    xlim(-200, 200)

    subplot(234)
    plot(gdt1, rad2deg(gphis_rec), ',', c='green')
    plot(dt1, rad2deg(phis_rec), ',', c='red')
    xlabel(r"$t_1 - t_3$ [ns]")
    ylabel(r"$\phi_{rec}$")
    xlim(-200, 200)

    subplot(235)
    plot(gdt2, rad2deg(gphis_rec), ',', c='green')
    plot(dt2, rad2deg(phis_rec), ',', c='red')
    xlabel(r"$t_1 - t_4$ [ns]")
    ylabel(r"$\phi_{rec}$")
    xlim(-200, 200)

    subplot(233)
    plot(gdt1, gdt2, ',', c='green')
    plot(dt1, dt2, ',', c='red')
    xlabel(r"$t_1 - t_3$ [ns]")
    ylabel(r"$t_1 - t_4$ [ns]")
    xlim(-200, 200)
    ylim(-200, 200)

    utils.saveplot()


def plot_failed_histograms():
    figure()

    global dt1, dt2, phis

    c = 3e8 * 1e-9
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    dt1 = array(dt1)
    dt2 = array(dt2)
    phis = array(phis)

    subplot(121)
    hist(c * dt1 / (10 * cos(phis - phi1)), bins=linspace(-20, 20, 100))
    xlabel(r"$c \, \Delta t_1 / (r_1 \cos(\phi - \phi_1))$")

    subplot(122)
    hist(c * dt2 / (10 * cos(phis - phi2)), bins=linspace(-20, 20, 100))
    xlabel(r"$c \, \Delta t_2 / (r_2 \cos(\phi - \phi_2))$")

    utils.saveplot()


def plot_uncertainty_zenith_angular_distance(group):
    group = group.E_1PeV
    rec = DirectionReconstruction

    N = 2

    # constants for uncertainty estimation
    # BEWARE: stations must be the same over all reconstruction tables used
    station = group.zenith_0.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    figure()
    graph = GraphArtist()
    # Uncertainty estimate
    x = linspace(0, deg2rad(45), 50)
    #x = array([pi / 8])
    phis = linspace(-pi, pi, 50)
    y, y2 = [], []
    for t in x:
        y.append(mean(rec.rel_phi_errorsq(t, phis, phi1, phi2, r1, r2)))
        y2.append(mean(rec.rel_theta1_errorsq(t, phis, phi1, phi2, r1, r2)))
    y = TIMING_ERROR * sqrt(array(y))
    y2 = TIMING_ERROR * sqrt(array(y2))
    ang_dist = sqrt((y * sin(x)) ** 2 + y2 ** 2)
    #plot(rad2deg(x), rad2deg(y), label="Estimate Phi")
    #plot(rad2deg(x), rad2deg(y2), label="Estimate Theta")
    plot(rad2deg(x), rad2deg(ang_dist), label="Angular distance")
    graph.plot(rad2deg(x), rad2deg(ang_dist), mark=None)
    print rad2deg(x)
    print rad2deg(y)
    print rad2deg(y2)
    print rad2deg(y * sin(x))
    print rad2deg(ang_dist)

    # Labels etc.
    xlabel("Shower zenith angle [deg]")
    ylabel("Angular distance [deg]")
    graph.set_xlabel(r"Shower zenith angle [\si{\degree}]")
    graph.set_ylabel(r"Angular distance [\si{\degree}]")
    graph.set_ylimits(min=6)
    #title(r"$N_{MIP} \geq %d$" % N)
    #ylim(0, 100)
    #legend(numpoints=1)
    utils.saveplot()
    artist.utils.save_graph(graph, dirname='plots')
    print


if __name__ == '__main__':
    # invalid values in arcsin will be ignored (nan handles the situation
    # quite well)
    np.seterr(invalid='ignore', divide='ignore')

    try:
        data
    except NameError:
        data = tables.open_file('master-ch4v2.h5', 'r')

    if '/reconstructions' not in data:
        print "Reconstructing shower direction..."
        do_full_reconstruction(data)
    else:
        print "Skipping reconstruction!"

    utils.set_prefix("DIR-")
    artist.utils.set_prefix("DIR-")
    do_reconstruction_plots(data)

    # These currently don't work
#    utils.set_prefix("WIP-")
#    do_jos_plots(data)
