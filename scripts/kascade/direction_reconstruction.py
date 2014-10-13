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

from myshowerfront import *

from artist import GraphArtist, MultiPlot
import artist.utils


DATADIR = '../simulations/plots'

USE_TEX = True

TIMING_ERROR = 2.4

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


def do_reconstruction_plots(data, table):
    """Make plots based upon earlier reconstructions"""

    plot_uncertainty_mip(table)
    plot_uncertainty_zenith(table)
    plot_uncertainty_core_distance(table)

    plot_phi_reconstruction_results_for_MIP(table, 1)
    plot_phi_reconstruction_results_for_MIP(table, 2)
    plot_theta_reconstruction_results_for_MIP(table, 1)
    plot_theta_reconstruction_results_for_MIP(table, 2)
    boxplot_theta_reconstruction_results_for_MIP(table, 1)
    boxplot_theta_reconstruction_results_for_MIP(table, 2)
    boxplot_phi_reconstruction_results_for_MIP(table, 1)
    boxplot_phi_reconstruction_results_for_MIP(table, 2)
    boxplot_arrival_times(table, 1)
    boxplot_core_distances_for_mips(table)


def do_lint_comparison(data):
    fsot = data.root.reconstructions_offsets
    lint = data.root.lint_reconstructions_offsets

    plot_fsot_vs_lint_for_zenith(fsot, lint)


def plot_uncertainty_mip(table):
    rec = DirectionReconstruction

    # constants for uncertainty estimation
    station = table.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    THETA = deg2rad(22.5)
    DTHETA = deg2rad(5.)
    DN = .1
    LOGENERGY = 15
    DLOGENERGY = .5

    figure()
    x, y, y2 = [], [], []
    for N in range(1, 6):
        x.append(N)
        events = table.read_where('(abs(min_n134 - N) <= DN) & (abs(reference_theta - THETA) <= DTHETA) & (abs(log10(k_energy) - LOGENERGY) <= DLOGENERGY)')
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
    print "mip: min_n134, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print

    # Simulation data
    sx, sy, sy2 = loadtxt(os.path.join(DATADIR, 'DIR-plot_uncertainty_mip.txt'))

    # Uncertainty estimate
    ex = linspace(1, 5, 50)
    phis = linspace(-pi, pi, 50)
    phi_errsq = mean(rec.rel_phi_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    theta_errsq = mean(rec.rel_theta1_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    #ey = TIMING_ERROR * std_t(ex) * sqrt(phi_errsq)
    #ey2 = TIMING_ERROR * std_t(ex) * sqrt(theta_errsq)

    R_list = [30, 20, 16, 14, 12]
    with tables.open_file('master-ch4v2.h5') as data2:
        mc = my_std_t_for_R(data2, x, R_list)
    mc = sqrt(mc ** 2 + 1.2 ** 2 + 2.5 ** 2)
    print mc
    ey = mc * sqrt(phi_errsq)
    ey2 = mc * sqrt(theta_errsq)

    nx = linspace(1, 4, 100)
    ey = spline(x, ey, nx)
    ey2 = spline(x, ey2, nx)

    # Plots
    plot(x, rad2deg(y), '^', label="Theta")
    plot(sx, rad2deg(sy), '^', label="Theta (sim)")
    plot(nx, rad2deg(ey2))#, label="Estimate Theta")
    plot(x, rad2deg(y2), 'v', label="Phi")
    plot(sx, rad2deg(sy2), 'v', label="Phi (sim)")
    plot(nx, rad2deg(ey))#, label="Estimate Phi")

    # Labels etc.
    xlabel("$N_{MIP} \pm %.1f$" % DN)
    ylabel("Angle reconstruction uncertainty [deg]")
    title(r"$\theta = 22.5^\circ \pm %d^\circ \quad %.1f \leq \log(E) \leq %.1f$" % (rad2deg(DTHETA), LOGENERGY - DLOGENERGY, LOGENERGY + DLOGENERGY))
    legend(numpoints=1)
    xlim(0.5, 4.5)
    utils.saveplot()
    print

    graph = GraphArtist()
    graph.plot(x, rad2deg(y), mark='o', linestyle=None)
    graph.plot(sx, rad2deg(sy), mark='square', linestyle=None)
    graph.plot(nx, rad2deg(ey2), mark=None)
    graph.plot(x, rad2deg(y2), mark='*', linestyle=None)
    graph.plot(sx, rad2deg(sy2), mark='square*', linestyle=None)
    graph.plot(nx, rad2deg(ey), mark=None)
    graph.set_xlabel(r"$N_\mathrm{MIP} \pm %.1f$" % DN)
    graph.set_ylabel(r"Angle reconstruction uncertainty [\si{\degree}]")
    graph.set_xlimits(max=4.5)
    graph.set_ylimits(0, 40)
    graph.set_xticks(range(5))
    artist.utils.save_graph(graph, dirname='plots')


def plot_uncertainty_zenith(table):
    rec = DirectionReconstruction

    # constants for uncertainty estimation
    station = table.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    N = 2
    DTHETA = deg2rad(1.)
    DN = .1
    LOGENERGY = 15
    DLOGENERGY = .5

    figure()
    rcParams['text.usetex'] = False
    x, y, y2 = [], [], []
    for theta in 5, 10, 15, 22.5, 30, 35:
        x.append(theta)
        THETA = deg2rad(theta)
        events = table.read_where('(min_n134 >= N) & (abs(reference_theta - THETA) <= DTHETA) & (abs(log10(k_energy) - LOGENERGY) <= DLOGENERGY)')
        print theta, len(events),
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
    print "zenith: theta, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print

    # Simulation data
    sx, sy, sy2 = loadtxt(os.path.join(DATADIR, 'DIR-plot_uncertainty_zenith.txt'))

    # Uncertainty estimate
    ex = linspace(0, deg2rad(35), 50)
    phis = linspace(-pi, pi, 50)
    ey, ey2, ey3 = [], [], []
    for t in ex:
        ey.append(mean(rec.rel_phi_errorsq(t, phis, phi1, phi2, r1, r2)))
        ey3.append(mean(rec.rel_phi_errorsq(t, phis, phi1, phi2, r1, r2)) * sin(t) ** 2)
        ey2.append(mean(rec.rel_theta1_errorsq(t, phis, phi1, phi2, r1, r2)))
    ey = TIMING_ERROR * sqrt(array(ey))
    ey3 = TIMING_ERROR * sqrt(array(ey3))
    ey2 = TIMING_ERROR * sqrt(array(ey2))

    graph = GraphArtist()

    # Plots
    plot(x, rad2deg(y), '^', label="Theta")
    graph.plot(x, rad2deg(y), mark='o', linestyle=None)
    #plot(sx, rad2deg(sy), '^', label="Theta (sim)")
    plot(rad2deg(ex), rad2deg(ey2))#, label="Estimate Theta")
    graph.plot(rad2deg(ex), rad2deg(ey2), mark=None)
    # Azimuthal angle undefined for zenith = 0
    plot(x[1:], rad2deg(y2[1:]), 'v', label="Phi")
    graph.plot(x[1:], rad2deg(y2[1:]), mark='*', linestyle=None)
    #plot(sx[1:], rad2deg(sy2[1:]), 'v', label="Phi (sim)")
    plot(rad2deg(ex), rad2deg(ey))#, label="Estimate Phi")
    graph.plot(rad2deg(ex), rad2deg(ey), mark=None)
    #plot(rad2deg(ex), rad2deg(ey3), label="Estimate Phi * sin(Theta)")

    # Labels etc.
    xlabel(r"Shower zenith angle [deg $\pm %d^\circ$]" % rad2deg(DTHETA))
    graph.set_xlabel(r"Shower zenith angle [\si{\degree}] $\pm \SI{%d}{\degree}$" % rad2deg(DTHETA))
    ylabel("Angle reconstruction uncertainty [deg]")
    graph.set_ylabel(r"Angle reconstruction uncertainty [\si{\degree}]")
    title(r"$N_{MIP} \geq %d, \quad %.1f \leq \log(E) \leq %.1f$" % (N, LOGENERGY - DLOGENERGY, LOGENERGY + DLOGENERGY))
    ylim(0, 60)
    graph.set_ylimits(0, 60)
    xlim(-.5, 37)
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    utils.saveplot()
    artist.utils.save_graph(graph, dirname='plots')
    print


def plot_uncertainty_core_distance(table):
    N = 2
    THETA = deg2rad(22.5)
    DTHETA = deg2rad(5.)
    DN = .5
    DR = 10
    LOGENERGY = 15
    DLOGENERGY = .5

    figure()
    x, y, y2 = [], [], []
    for R in range(0, 81, 20):
        x.append(R)
        events = table.read_where('(abs(min_n134 - N) <= DN) & (abs(reference_theta - THETA) <= DTHETA) & (abs(r - R) <= DR) & (abs(log10(k_energy) - LOGENERGY) <= DLOGENERGY)')
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

#    # Simulation data
    sx, sy, sy2 = loadtxt(os.path.join(DATADIR, 'DIR-plot_uncertainty_core_distance.txt'))

    graph = GraphArtist()

    # Plots
    plot(x, rad2deg(y), '^-', label="Theta")
    graph.plot(x[:-1], rad2deg(y[:-1]), mark='o')
    plot(sx, rad2deg(sy), '^-', label="Theta (sim)")
    graph.plot(sx[:-1], rad2deg(sy[:-1]), mark='square')
    plot(x, rad2deg(y2), 'v-', label="Phi")
    graph.plot(x[:-1], rad2deg(y2[:-1]), mark='*')
    plot(sx, rad2deg(sy2), 'v-', label="Phi (sim)")
    graph.plot(sx[:-1], rad2deg(sy2[:-1]), mark='square*')

    # Labels etc.
    xlabel("Core distance [m] $\pm %d$" % DR)
    graph.set_xlabel(r"Core distance [\si{\meter}] $\pm \SI{%d}{\meter}$" % DR)
    ylabel("Angle reconstruction uncertainty [deg]")
    graph.set_ylabel(r"Angle reconstruction uncertainty [\si{\degree}]")
    title(r"$N_{MIP} = %d \pm %.1f, \theta = 22.5^\circ \pm %d^\circ, %.1f \leq \log(E) \leq %.1f$" % (N, DN, rad2deg(DTHETA), LOGENERGY - DLOGENERGY, LOGENERGY + DLOGENERGY))
    ylim(ymin=0)
    graph.set_ylimits(min=0)
    xlim(-2, 62)
    legend(numpoints=1, loc='best')
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


def plot_phi_reconstruction_results_for_MIP(table, N):
    THETA = deg2rad(22.5)
    DTHETA = deg2rad(5.)

    events = table.read_where('(min_n134 >= N) & (abs(reference_theta - THETA) <= DTHETA)')
    sim_phi = events['reference_phi']
    r_phi = events['reconstructed_phi']

    figure()
    plot_2d_histogram(rad2deg(sim_phi), rad2deg(r_phi), 180)
    xlabel(r"$\phi_K$ [deg]")
    ylabel(r"$\phi_H$ [deg]")
    title(r"$N_{MIP} \geq %d, \quad \theta = 22.5^\circ \pm %d^\circ$" % (N, rad2deg(DTHETA)))

    utils.saveplot(N)

    graph = artist.GraphArtist()
    bins = linspace(-180, 180, 73)
    H, x_edges, y_edges = histogram2d(rad2deg(sim_phi), rad2deg(r_phi),
                                      bins=bins)
    graph.histogram2d(H, x_edges, y_edges, type='reverse_bw')
    graph.set_xlabel(r'$\phi_K$ [\si{\degree}]')
    graph.set_ylabel(r'$\phi_H$ [\si{\degree}]')
    graph.set_xticks(range(-180, 181, 90))
    graph.set_yticks(range(-180, 181, 90))
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def plot_theta_reconstruction_results_for_MIP(table, N):
    events = table.read_where('min_n134 >= N')
    sim_theta = events['reference_theta']
    r_theta = events['reconstructed_theta']

    figure()
    x_edges = linspace(0, 40, 81)
    y_edges = linspace(0, 40, 81)
    plot_2d_histogram(rad2deg(sim_theta), rad2deg(r_theta), (x_edges, y_edges))
    xlabel(r"$\theta_K$ [deg]")
    ylabel(r"$\theta_H$ [deg]")
    title(r"$N_{MIP} \geq %d$" % N)

    utils.saveplot(N)

    graph = artist.GraphArtist()
    bins = linspace(0, 40, 41)
    H, x_edges, y_edges = histogram2d(rad2deg(sim_theta), rad2deg(r_theta),
                                      bins=bins)
    graph.histogram2d(H, x_edges, y_edges, type='reverse_bw')
    graph.set_xlabel(r'$\theta_K$ [\si{\degree}]')
    graph.set_ylabel(r'$\theta_H$ [\si{\degree}]')
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def boxplot_theta_reconstruction_results_for_MIP(table, N):
    figure()

    DTHETA = deg2rad(1.)

    angles = [0, 5, 10, 15, 22.5, 35]
    r_dtheta = []
    x = []
    d25, d50, d75 = [], [], []
    for angle in angles:
        theta = deg2rad(angle)
        sel = table.read_where('(min_n134 >= N) & (abs(reference_theta - theta) <= DTHETA)')
        dtheta = rad2deg(sel[:]['reconstructed_theta'] - sel[:]['reference_theta'])
        r_dtheta.append(dtheta)

        d25.append(scoreatpercentile(dtheta, 25))
        d50.append(scoreatpercentile(dtheta, 50))
        d75.append(scoreatpercentile(dtheta, 75))
        x.append(angle)

    #boxplot(r_dtheta, sym='', positions=angles, widths=2.)
    fill_between(x, d25, d75, color='0.75')
    plot(x, d50, 'o-', color='black')

    xlabel(r"$\theta_K$ [deg]")
    ylabel(r"$\theta_H - \theta_K$ [deg]")
    title(r"$N_{MIP} \geq %d$" % N)

    axhline(0, color='black')
    ylim(-20, 25)
    xlim(0, 35)

    utils.saveplot(N)

    graph = GraphArtist()
    graph.draw_horizontal_line(0, linestyle='gray')
    graph.shade_region(angles, d25, d75)
    graph.plot(angles, d50, linestyle=None)
    graph.set_xlabel(r"$\theta_K$ [\si{\degree}]")
    graph.set_ylabel(r"$\theta_H - \theta_K$ [\si{\degree}]")
    graph.set_ylimits(-5, 15)
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def boxplot_phi_reconstruction_results_for_MIP(table, N):
    figure()

    THETA = deg2rad(22.5)
    DTHETA = deg2rad(5.)

    bin_edges = linspace(-180, 180, 18)
    x, r_dphi = [], []
    d25, d50, d75 = [], [], []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        rad_low = deg2rad(low)
        rad_high = deg2rad(high)
        query = '(min_n134 >= N) & (rad_low < reference_phi) & (reference_phi < rad_high) & (abs(reference_theta - THETA) <= DTHETA)'
        sel = table.read_where(query)
        dphi = sel[:]['reconstructed_phi'] - sel[:]['reference_phi']
        dphi = (dphi + pi) % (2 * pi) - pi
        r_dphi.append(rad2deg(dphi))

        d25.append(scoreatpercentile(rad2deg(dphi), 25))
        d50.append(scoreatpercentile(rad2deg(dphi), 50))
        d75.append(scoreatpercentile(rad2deg(dphi), 75))
        x.append((low + high) / 2)

    #boxplot(r_dphi, positions=x, widths=1 * (high - low), sym='')
    fill_between(x, d25, d75, color='0.75')
    plot(x, d50, 'o-', color='black')

    xlabel(r"$\phi_K$ [deg]")
    ylabel(r"$\phi_H - \phi_K$ [deg]")
    title(r"$N_{MIP} \geq %d, \quad \theta = 22.5^\circ \pm %d^\circ$" % (N, rad2deg(DTHETA)))

    xticks(linspace(-180, 180, 9))
    axhline(0, color='black')

    utils.saveplot(N)

    graph = GraphArtist()
    graph.draw_horizontal_line(0, linestyle='gray')
    graph.shade_region(x, d25, d75)
    graph.plot(x, d50, linestyle=None)
    graph.set_xlabel(r"$\phi_K$ [\si{\degree}]")
    graph.set_ylabel(r"$\phi_H - \phi_K$ [\si{\degree}]")
    graph.set_xticks([-180, -90, '...', 180])
    graph.set_xlimits(-180, 180)
    graph.set_ylimits(-23, 23)
    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def boxplot_arrival_times(table, N):
    THETA = deg2rad(0)
    DTHETA = deg2rad(10.)

    LOGENERGY = 15
    DLOGENERGY = .5

    bin_edges = linspace(0, 80, 5)
    x = []
    t25, t50, t75 = [], [], []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        query = '(min_n134 >= N) & (low <= r) & (r < high) & (abs(reference_theta - THETA) <= DTHETA) & (abs(log10(k_energy) - LOGENERGY) <= DLOGENERGY)'
        sel = table.read_where(query)
        t2 = sel[:]['t2']
        t1 = sel[:]['t1']
        ct1 = t1.compress((t1 > -999) & (t2 > -999))
        ct2 = t2.compress((t1 > -999) & (t2 > -999))
        print len(t1), len(t2), len(ct1), len(ct2)
        t25.append(scoreatpercentile(abs(ct2 - ct1), 25))
        t50.append(scoreatpercentile(abs(ct2 - ct1), 50))
        t75.append(scoreatpercentile(abs(ct2 - ct1), 75))
        x.append((low + high) / 2)

    sx, st25, st50, st75 = loadtxt(os.path.join(DATADIR, 'DIR-boxplot_arrival_times-1.txt'))

    fig = figure()

    ax1 = subplot(131)
    fill_between(sx, st25, st75, color='0.75')
    plot(sx, st50, 'o-', color='black', markerfacecolor='none')

    ax2 = subplot(132, sharex=ax1, sharey=ax1)
    fill_between(x, t25, t75, color='0.75')
    plot(x, t50, 'o-', color='black')

    ax3 = subplot(133, sharex=ax1, sharey=ax1)
    plot(sx, st50, 'o-', color='black', markerfacecolor='none')
    plot(x, t50, 'o-', color='black')

    ax2.xaxis.set_label_text("Core distance [m]")
    ax1.yaxis.set_label_text("Arrival time difference $|t_2 - t_1|$ [ns]")
    fig.suptitle(r"$N_{MIP} \geq %d, \quad \theta = 0^\circ \pm %d^\circ, \quad %.1f \leq \log(E) \leq %.1f$" % (N, rad2deg(DTHETA), LOGENERGY - DLOGENERGY, LOGENERGY + DLOGENERGY))

    ylim(ymax=15)
    xlim(xmax=80)
    locator_params(tight=True, nbins=4)

    fig.subplots_adjust(left=.1, right=.95)

    fig.set_size_inches(5, 2.67)
    utils.saveplot(N)

    sx = sx.compress(sx < 80)
    st25 = st25[:len(sx)]
    st50 = st50[:len(sx)]
    st75 = st75[:len(sx)]

    graph = MultiPlot(1, 3, width=r'.3\linewidth', height=r'.4\linewidth')

    graph.shade_region(0, 0, sx, st25, st75)
    graph.plot(0, 0, sx, st50, linestyle=None)
    graph.plot(0, 2, sx, st50)
    graph.set_label(0, 0, 'simulation', 'upper left')

    graph.shade_region(0, 1, x, t25, t75)
    graph.plot(0, 1, x, t50, linestyle=None, mark='*')
    graph.plot(0, 2, x, t50, mark='*')
    graph.set_label(0, 1, 'experiment', 'upper left')

    graph.set_label(0, 2, 'sim + exp', 'upper left')

    graph.set_ylimits(0, 15)
    graph.set_xlimits(0, 80)
    graph.set_xlabel("Core distance [\si{\meter}]")
    graph.set_ylabel(r"Arrival time difference $|t_2 - t_1|$ [\si{\nano\second}]")
    graph.show_xticklabels_for_all([(0, 0), (0, 1), (0, 2)])
    graph.set_xticklabels_position(0, 1, 'right')
    graph.show_yticklabels(0, 0)

    artist.utils.save_graph(graph, suffix=N, dirname='plots')


def boxplot_core_distances_for_mips(table):
    THETA = deg2rad(22.5)
    DTHETA = deg2rad(1.)
    DN = .5

    ENERGY = 1e15
    DENERGY = 2e14

    MAX_R = 80

    r25_list = []
    r50_list = []
    r75_list = []
    x = []
    for N in range(1, 5):
        sel = table.read_where('(abs(min_n134 - N) <= DN) & (abs(reference_theta - THETA) <= DTHETA) & (abs(k_energy - ENERGY) <= DENERGY) & (r <= MAX_R)')
        r = sel[:]['r']
        r25_list.append(scoreatpercentile(r, 25))
        r50_list.append(scoreatpercentile(r, 50))
        r75_list.append(scoreatpercentile(r, 75))
        x.append(N)
        print len(r)

    sx, sr25, sr50, sr75 = loadtxt(os.path.join(DATADIR, 'DIR-save_for_kascade_boxplot_core_distances_for_mips.txt'))

    fig = figure()

    ax1 = subplot(131)
    fill_between(sx, sr25, sr75, color='0.75')
    plot(sx, sr50, 'o-', color='black', markerfacecolor='none')

    ax2 = subplot(132, sharex=ax1, sharey=ax1)
    fill_between(x, r25_list, r75_list, color='0.75')
    plot(x, r50_list, 'o-', color='black')

    ax3 = subplot(133, sharex=ax1, sharey=ax1)
    plot(x, sr50, 'o-', color='black', markerfacecolor='none')
    plot(x, r50_list, 'o-', color='black')

    ax2.xaxis.set_label_text("Minimum number of particles $\pm %.1f$" % DN)
    ax1.yaxis.set_label_text("Core distance [m]")
    fig.suptitle(r"$\theta = 22.5^\circ \pm %d^\circ, \quad %.1f \leq \log(E) \leq %.1f$" % (rad2deg(DTHETA), log10(ENERGY - DENERGY), log10(ENERGY + DENERGY)))

    ax1.xaxis.set_ticks([1, 2, 3, 4])
    fig.subplots_adjust(left=.1, right=.95)

    ylim(ymin=0)
    xlim(.8, 4.2)

    fig.set_size_inches(5, 2.67)
    utils.saveplot()

    graph = MultiPlot(1, 3, width=r'.3\linewidth', height=r'.4\linewidth')

    graph.shade_region(0, 0, sx, sr25, sr75)
    graph.plot(0, 0, sx, sr50, linestyle=None)
    graph.plot(0, 2, sx, sr50)
    graph.set_label(0, 0, 'simulation')

    graph.shade_region(0, 1, x, r25_list, r75_list)
    graph.plot(0, 1, x, r50_list, linestyle=None, mark='*')
    graph.plot(0, 2, x, r50_list, mark='*')
    graph.set_label(0, 1, 'experiment')

    graph.set_label(0, 2, 'sim + exp')

    graph.set_ylimits(0, 50)
    graph.set_xlabel("Minimum number of particles $\pm %.1f$" % DN)
    graph.set_ylabel("Core distance [\si{\meter}]")
    graph.show_xticklabels_for_all([(0, 0), (0, 1), (0, 2)])
    graph.show_yticklabels(0, 0)

    artist.utils.save_graph(graph, dirname='plots')


def plot_2d_histogram(x, y, bins):
    H, xedges, yedges = histogram2d(x, y, bins)
    imshow(H.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
           origin='lower left', interpolation='lanczos', aspect='auto',
           cmap=cm.Greys)
    colorbar()


def plot_fsot_vs_lint_for_zenith(fsot, lint):
    bins = linspace(0, 35, 21)

    min_N = 1

    x, f_y, f_y2, l_y, l_y2 = [], [], [], [], []
    for low, high in zip(bins[:-1], bins[1:]):
        rad_low = deg2rad(low)
        rad_high = deg2rad(high)

        query = '(min_n134 >= min_N) & (rad_low <= reference_theta) & (reference_theta < rad_high)'
        f_sel = fsot.read_where(query)
        l_sel = lint.read_where(query)

        errors = f_sel['reconstructed_phi'] - f_sel['reference_phi']
        errors2 = f_sel['reconstructed_theta'] - f_sel['reference_theta']
        #f_y.append(std(errors))
        #f_y2.append(std(errors2))
        f_y.append((scoreatpercentile(errors, 83) - scoreatpercentile(errors, 17)) / 2)
        f_y2.append((scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17)) / 2)

        errors = l_sel['reconstructed_phi'] - l_sel['reference_phi']
        errors2 = l_sel['reconstructed_theta'] - l_sel['reference_theta']
        #l_y.append(std(errors))
        #l_y2.append(std(errors2))
        l_y.append((scoreatpercentile(errors, 83) - scoreatpercentile(errors, 17)) / 2)
        l_y2.append((scoreatpercentile(errors2, 83) - scoreatpercentile(errors2, 17)) / 2)

        x.append((low + high) / 2)

        print x[-1], len(f_sel), len(l_sel)

    clf()
    plot(x, rad2deg(f_y), label="FSOT phi")
    plot(x, rad2deg(f_y2), label="FSOT theta")
    plot(x, rad2deg(l_y), label="LINT phi")
    plot(x, rad2deg(l_y2), label="LINT theta")
    legend()
    xlabel("Shower zenith angle [deg]")
    ylabel("Angle reconstruction uncertainty [deg]")
    title(r"$N_{MIP} \geq %d$" % min_N)
    utils.saveplot()

    graph = GraphArtist()
    graph.plot(x, rad2deg(f_y), mark=None)
    graph.plot(x, rad2deg(l_y), mark=None, linestyle='dashed')
    graph.plot(x, rad2deg(f_y2), mark=None)
    graph.plot(x, rad2deg(l_y2), mark=None, linestyle='dashed')
    graph.set_xlabel(r"Shower zenith angle [\si{\degree}]")
    graph.set_ylabel(r"Angle reconstruction uncertainty [\si{\degree}]")
    artist.utils.save_graph(graph, dirname='plots')


if __name__ == '__main__':
    # invalid values in arcsin will be ignored (nan handles the situation
    # quite well)
    np.seterr(invalid='ignore', divide='ignore')

    try:
        data
    except NameError:
        data = tables.open_file('kascade.h5', 'r')

    artist.utils.set_prefix("KAS-")
    utils.set_prefix("KAS-")
    do_reconstruction_plots(data, data.root.reconstructions)
    do_lint_comparison(data)
    artist.utils.set_prefix("KAS-LINT-")
    utils.set_prefix("KAS-LINT-")
    do_reconstruction_plots(data, data.root.lint_reconstructions)
    artist.utils.set_prefix("KAS-OFFSETS-")
    utils.set_prefix("KAS-OFFSETS-")
    do_reconstruction_plots(data, data.root.reconstructions_offsets)
    artist.utils.set_prefix("KAS-LINT-OFFSETS-")
    utils.set_prefix("KAS-LINT-OFFSETS-")
    do_reconstruction_plots(data, data.root.lint_reconstructions_offsets)
