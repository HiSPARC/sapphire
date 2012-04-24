from __future__ import division

import tables
from itertools import izip
import os.path

import progressbar as pb

from pylab import *

from scipy import integrate
from scipy.special import erf
from scipy.stats import scoreatpercentile

import utils

from sapphire.analysis import DirectionReconstruction, BinnedDirectionReconstruction


USE_TEX = True

TIMING_ERROR = 2.8

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
        data.removeNode('/reconstructions', recursive=True)

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
    plot_reconstruction_efficiency_vs_R_for_angles(1)
    plot_reconstruction_efficiency_vs_R_for_angles(2)
    plot_reconstruction_efficiency_vs_R_for_mips()

def plot_uncertainty_mip(group):
    table = group.E_1PeV.zenith_22_5
    rec = DirectionReconstruction

    # constants for uncertainty estimation
    station = table.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    figure()
    x, y, y2 = [], [], []
    for N in range(1, 6):
        x.append(N)
        events = table.readWhere('min_n134==%d' % N)
        print len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))

    plot(x, rad2deg(y), '^', label="Theta")
    plot(x, rad2deg(y2), 'v', label="Phi")
    print
    print "mip: min_n134, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    utils.savedata((x, y, y2))

    # Uncertainty estimate
    x = linspace(1, 5, 50)
    phis = linspace(-pi, pi, 50)
    phi_errsq = mean(rec.rel_phi_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    theta_errsq = mean(rec.rel_theta1_errorsq(pi / 8, phis, phi1, phi2, r1, r2))
    y = TIMING_ERROR * std_t(x) * sqrt(phi_errsq)
    y2 = TIMING_ERROR * std_t(x) * sqrt(theta_errsq)
    plot(x, rad2deg(y), label="Estimate Phi")
    plot(x, rad2deg(y2), label="Estimate Theta")
    # Labels etc.
    xlabel("Minimum number of particles")
    ylabel("Angle reconstruction uncertainty [deg]")
    #title(r"$\theta = 22.5^\circ$")
    legend(numpoints=1)
    utils.saveplot()
    print

def plot_uncertainty_zenith(group):
    group = group.E_1PeV
    rec = DirectionReconstruction

    N = 2

    # constants for uncertainty estimation
    # BEWARE: stations must be the same over all reconstruction tables used
    station = group.zenith_0.attrs.cluster.stations[0]
    r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
    r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

    figure()
    x, y, y2 = [], [], []
    for THETA in 0, 5, 10, 15, 22.5, 30, 35, 45:
        x.append(THETA)
        table = group._f_getChild('zenith_%s' % str(THETA).replace('.', '_'))
        events = table.readWhere('min_n134 >= N')
        print THETA, len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(x, rad2deg(y), '^', label="Theta")
    # Azimuthal angle undefined for zenith = 0
    plot(x[1:], rad2deg(y2[1:]), 'v', label="Phi")
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
    plot(rad2deg(x), rad2deg(y2), label="Estimate Theta")

    # Labels etc.
    xlabel("Shower zenith angle [deg]")
    ylabel("Angle reconstruction uncertainty [deg]")
    #title(r"$N_{MIP} \geq %d$" % N)
    ylim(0, 100)
    legend(numpoints=1)
    utils.saveplot()
    print

def plot_uncertainty_core_distance(group):
    table = group.E_1PeV.zenith_22_5

    N = 2
    DR = 10

    figure()
    x, y, y2 = [], [], []
    for R in range(0, 81, 20):
        x.append(R)
        events = table.readWhere('(min_n134 == N) & (abs(r - R) <= DR)')
        print len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))

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
            table = group._f_getChild('zenith_22_5_size%d' % size)
        else:
            table = group._f_getChild('zenith_22_5')

        events = table.readWhere('min_n134 >= N')
        print size, len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(x, rad2deg(y), '^', label="Theta")
    plot(x, rad2deg(y2), 'v', label="Phi")
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
    plot(x, rad2deg(y2), label="Estimate Theta")

    # Labels etc.
    xlabel("Station size [m]")
    ylabel("Angle reconstruction uncertainty [deg]")
    #title(r"$\theta = 22.5^\circ, N_{MIP} \geq %d$" % N)
    legend(numpoints=1)
    utils.saveplot()
    print

def plot_uncertainty_binsize(group):
    group = group.E_1PeV
    rec = DirectionReconstruction

    N = 2

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
            table = group._f_getChild('zenith_22_5_binned_randomized_%s' % str(bin_size).replace('.', '_'))
        else:
            table = group.zenith_22_5
        events = table.readWhere('min_n134 >= 2')

        print bin_size, len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(x, rad2deg(y), '^', label="Theta")
    plot(x, rad2deg(y2), 'v', label="Phi")
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
    plot(x, rad2deg(y2), label="Estimate Theta")

    # Labels etc.
    xlabel("Bin size [ns]")
    ylabel("Angle reconstruction uncertainty [deg]")
    #title(r"$\theta = 22.5^\circ, N_{MIP} \geq %d$" % N)
    legend(loc='best', numpoints=1)
    ylim(ymin=0)
    utils.saveplot()
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

    events = table.readWhere('min_n134 >= %d' % N)
    sim_phi = events['reference_phi']
    r_phi = events['reconstructed_phi']

    figure()
    plot_2d_histogram(rad2deg(sim_phi), rad2deg(r_phi), 180)
    xlabel(r"$\phi_{simulated}$ [deg]")
    ylabel(r"$\phi_{reconstructed}$ [deg]")
    #title(r"$N_{MIP} \geq %d, \quad \theta = 22.5^\circ$" % N)

    utils.saveplot(N)

def boxplot_theta_reconstruction_results_for_MIP(group, N):
    group = group.E_1PeV

    figure()

    angles = [0, 5, 10, 15, 22.5, 30, 35, 45]
    r_dtheta = []
    d25, d50, d75 = [], [], []
    for angle in angles:
        table = group._f_getChild('zenith_%s' % str(angle).replace('.', '_'))
        sel = table.readWhere('min_n134 >= %d' % N)
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
        sel = table.readWhere(query)
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

def boxplot_arrival_times(group, N):
    table = group.E_1PeV.zenith_0

    sel = table.readWhere('min_n134 >= N')
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
        sel = table.readWhere(query)
        t1 = sel[:]['t1']
        t3 = sel[:]['t3']
        t4 = sel[:]['t4']
        ts = concatenate([t1, t3, t4])

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

def boxplot_core_distances_for_mips(group):
    table = group.E_1PeV.zenith_22_5
    
    figure()

    r_list = []
    r25, r50, r75 = [], [], []
    x = []
    for N in range(1, 5):
        sel = table.readWhere('min_n134 == N')
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

def save_for_kascade_boxplot_core_distances_for_mips(group):
    table = group.E_1PeV.zenith_22_5

    r25_list = []
    r50_list = []
    r75_list = []
    x = []
    for N in range(1, 5):
        sel = table.readWhere('(min_n134 == N) & (r <= 80)')
        r = sel[:]['r']
        r25_list.append(scoreatpercentile(r, 25))
        r50_list.append(scoreatpercentile(r, 50))
        r75_list.append(scoreatpercentile(r, 75))
        x.append(N)

    utils.savedata((x, r25_list, r50_list, r75_list))

def plot_detection_efficiency_vs_R_for_angles(N):
    figure()

    bin_edges = linspace(0, 100, 20)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2.

    for angle in [0, 22.5, 35]:
        angle_str = str(angle).replace('.', '_')
        shower_group = '/simulations/E_1PeV/zenith_%s' % angle_str

        efficiencies = []
        for low, high in zip(bin_edges[:-1], bin_edges[1:]):
            shower_results = []
            for shower in data.listNodes(shower_group):
                sel_query = '(low <= r) & (r < high)'
                coinc_sel = shower.coincidences.readWhere(sel_query)
                ids = coinc_sel['id']
                obs_sel = shower.observables.readCoordinates(ids)
                assert (obs_sel['id'] == ids).all()

                o = obs_sel
                sel = obs_sel.compress((o['n1'] >= N) & (o['n3'] >= N) &
                                       (o['n4'] >= N))
                shower_results.append(len(sel) / len(obs_sel))
            efficiencies.append(mean(shower_results))

        plot(x, efficiencies, label=r'$\theta = %s^\circ$' % angle)

    xlabel("Core distance [m]")
    ylabel("Detection efficiency")
    #title(r"$N_{MIP} \geq %d$" % N)
    legend()

    utils.saveplot(N)

def plot_reconstruction_efficiency_vs_R_for_angles(N):
    group = data.root.reconstructions.E_1PeV

    figure()

    bin_edges = linspace(0, 100, 10)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2.

    for angle in [0, 22.5, 35]:
        angle_str = str(angle).replace('.', '_')
        shower_group = '/simulations/E_1PeV/zenith_%s' % angle_str
        reconstructions = group._f_getChild('zenith_%s' % angle_str)

        efficiencies = []
        for low, high in zip(bin_edges[:-1], bin_edges[1:]):
            shower_results = []
            for shower in data.listNodes(shower_group):
                sel_query = '(low <= r) & (r < high)'
                coinc_sel = shower.coincidences.readWhere(sel_query)
                ids = coinc_sel['id']
                obs_sel = shower.observables.readCoordinates(ids)
                assert (obs_sel['id'] == ids).all()

                o = obs_sel
                sel = obs_sel.compress((o['n1'] >= N) & (o['n3'] >= N) &
                                       (o['n4'] >= N))
                shower_results.append(len(sel))
            ssel = reconstructions.readWhere('(min_n134 >= N) & (low <= r) & (r < high)')
            efficiencies.append(len(ssel) / sum(shower_results))

        plot(x, efficiencies, label=r'$\theta = %s^\circ$' % angle)

    xlabel("Core distance [m]")
    ylabel("Reconstruction efficiency")
    #title(r"$N_{MIP} \geq %d$" % N)
    legend()

    utils.saveplot(N)

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
            for shower in data.listNodes(shower_group):
                sel_query = '(low <= r) & (r < high)'
                coinc_sel = shower.coincidences.readWhere(sel_query)
                ids = coinc_sel['id']
                obs_sel = shower.observables.readCoordinates(ids)
                assert (obs_sel['id'] == ids).all()

                o = obs_sel
                sel = o.compress(amin(array([o['n1'], o['n3'], o['n4']]), 0) == N)

                shower_results.append(len(sel))
            ssel = reconstructions.readWhere('(min_n134 == N) & (low <= r) & (r < high)')
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
           origin='lower left', interpolation='lanczos', aspect='auto')
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

    progressbar = pb.ProgressBar(maxval=len(observables))
    dt1, dt2, phis_sim, phis_rec = [], [], [], []
    gdt1, gdt2, gphis_sim, gphis_rec = [], [], [], []

    for event, coincidence in progressbar(izip(observables, coincidences)):
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


if __name__ == '__main__':
    # invalid values in arcsin will be ignored (nan handles the situation
    # quite well)
    np.seterr(invalid='ignore', divide='ignore')

    try:
        data
    except NameError:
        data = tables.openFile('master-ch4v2.h5', 'r')

    if '/reconstructions' not in data:
        print "Reconstructing shower direction..."
        do_full_reconstruction(data)
    else:
        print "Skipping reconstruction!"

    utils.set_prefix("DIR-")
    do_reconstruction_plots(data)

    # These currently don't work
#    utils.set_prefix("WIP-")
#    do_jos_plots(data)
