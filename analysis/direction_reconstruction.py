from __future__ import division

import tables
from itertools import combinations, izip
import re
import csv

import progressbar as pb

from pylab import *
from scipy.optimize import curve_fit
#from tikz_plot import tikz_2dhist

from scipy import integrate
from scipy.special import erf

import utils


USE_TEX = False

TIMING_ERROR = 4
#TIMING_ERROR = 7

# For matplotlib plots
if USE_TEX:
    rcParams['font.serif'] = 'Computer Modern'
    rcParams['font.sans-serif'] = 'Computer Modern'
    rcParams['font.family'] = 'sans-serif'
    rcParams['figure.figsize'] = [5 * x for x in (1, 2. / 3)]
    rcParams['figure.subplot.left'] = 0.125
    rcParams['figure.subplot.bottom'] = 0.125
    rcParams['font.size'] = 11
    rcParams['legend.fontsize'] = 'small'


class ReconstructedEvent(tables.IsDescription):
    r = tables.Float32Col()
    phi = tables.Float32Col()
    alpha = tables.Float32Col()
    t1 = tables.Float32Col()
    t2 = tables.Float32Col()
    t3 = tables.Float32Col()
    t4 = tables.Float32Col()
    n1 = tables.UInt16Col()
    n2 = tables.UInt16Col()
    n3 = tables.UInt16Col()
    n4 = tables.UInt16Col()
    reference_theta = tables.Float32Col()
    reference_phi = tables.Float32Col()
    reconstructed_theta = tables.Float32Col()
    reconstructed_phi = tables.Float32Col()
    min_n134 = tables.UInt16Col()
    size = tables.UInt8Col()
    bin = tables.Float32Col()
    bin_r = tables.BoolCol()


class DirectionReconstruction():
    def __init__(self, datafile, results_table, min_n134=1., N=None):
        self.data = datafile
        self.results_table = results_table
        self.min_n134 = min_n134
        self.N = N

    def reconstruct_angles(self, tablename, THETA, binning=False,
                           randomize_binning=False):
        """Reconstruct angles from simulation for minimum particle density"""

        shower_group = self.data.getNode('/simulations/E_1PeV', tablename)

        progressbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                              pb.ETA()],
                                     fd=sys.stderr)

        for shower in progressbar(self.data.listNodes(shower_group)):
            shower_table = shower.observables
            coincidence_table = shower.coincidences
            self.station, = shower._v_attrs.cluster.stations
            dst_row = self.results_table.row

            for event, coincidence in zip(shower_table[:self.N], coincidence_table[:self.N]):
                assert event['id'] == coincidence['id']
                if min(event['n1'], event['n3'], event['n4']) >= self.min_n134:
                    # Do we need to bin timing data?
                    if binning is not False:
                        event['t1'] = floor(event['t1'] / binning) * binning
                        event['t2'] = floor(event['t2'] / binning) * binning
                        event['t3'] = floor(event['t3'] / binning) * binning
                        event['t4'] = floor(event['t4'] / binning) * binning
                        # Do we need to randomize inside a bin?
                        if randomize_binning is True:
                            event['t1'] += uniform(0, binning)
                            event['t2'] += uniform(0, binning)
                            event['t3'] += uniform(0, binning)
                            event['t4'] += uniform(0, binning)

                    theta, phi = self.reconstruct_angle(event)

                    alpha = event['alpha']

                    if not isnan(theta) and not isnan(phi):
    #                    ang_dist = arccos(sin(theta) * sin(THETA) *
    #                                      cos(phi - -alpha) + cos(theta) *
    #                                      cos(THETA))

                        dst_row['r'] = coincidence['r']
                        dst_row['phi'] = coincidence['phi']
                        dst_row['alpha'] = event['alpha']
                        dst_row['t1'] = event['t1']
                        dst_row['t2'] = event['t2']
                        dst_row['t3'] = event['t3']
                        dst_row['t4'] = event['t4']
                        dst_row['n1'] = event['n1']
                        dst_row['n2'] = event['n2']
                        dst_row['n3'] = event['n3']
                        dst_row['n4'] = event['n4']
                        dst_row['reference_theta'] = THETA
                        dst_row['reference_phi'] = coincidence['alpha']
                        dst_row['reconstructed_theta'] = theta
                        dst_row['reconstructed_phi'] = phi
                        dst_row['min_n134'] = min(event['n1'], event['n3'], event['n4'])
                        r, phi = self.calc_r_and_phi(1, 3)
                        dst_row['size'] = r
                        if binning is False:
                            bin_size = 0
                        else:
                            bin_size = binning
                        dst_row['bin'] = bin_size
                        dst_row['bin_r'] = randomize_binning
                        dst_row.append()
            self.results_table.flush()

    def reconstruct_angle(self, event):
        """Reconstruct angles from a single event"""

        dt1 = event['t1'] - event['t3']
        dt2 = event['t1'] - event['t4']

        return self.reconstruct_angle_dt(dt1, dt2)

    def reconstruct_angle_dt(self, dt1, dt2):
        """Reconstruct angle given time differences"""

        c = 3.00e+8

        r1, phi1 = self.calc_r_and_phi(1, 3)
        r2, phi2 = self.calc_r_and_phi(1, 4)

        phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                      (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
        theta1 = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
        theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))

        e1 = sqrt(self.rel_theta1_errorsq(theta1, phi, phi1, phi2, r1, r2))
        e2 = sqrt(self.rel_theta2_errorsq(theta2, phi, phi1, phi2, r1, r2))

        theta_wgt = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

        return theta_wgt, phi

    def calc_r_and_phi(self, s1, s2):
        """Calculate angle between detectors (phi1, phi2)"""

        x1, y1 = self.station.detectors[s1 - 1].get_xy_coordinates()
        x2, y2 = self.station.detectors[s2 - 1].get_xy_coordinates()

        r = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        phi = arctan2((y2 - y1), (x2 - x1))

        return r, phi

    def rel_theta1_errorsq(self, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        sintheta = sin(theta)
        sinphiphi1 = sin(phi - phi1)

        den = (1 - sintheta ** 2) * r1 ** 2 * cos(phi - phi1) ** 2

        A = (r1 ** 2 * sinphiphi1 ** 2
             * self.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        B = (r1 * c * sinphiphi1
             * (self.dphi_dt0(theta, phi, phi1, phi2, r1, r2)
                - self.dphi_dt1(theta, phi, phi1, phi2, r1, r2)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 - 2 * B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    def rel_theta2_errorsq(self, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        sintheta = sin(theta)
        sinphiphi2 = sin(phi - phi2)

        den = (1 - sintheta ** 2) * r2 ** 2 * cos(phi - phi2) ** 2

        A = (r2 ** 2 * sinphiphi2 ** 2
             * self.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        B = (r2 * c * sinphiphi2
             * (self.dphi_dt0(theta, phi, phi1, phi2, r1, r2)
                - self.dphi_dt2(theta, phi, phi1, phi2, r1, r2)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 - 2 * B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    def rel_phi_errorsq(self, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) ** 2 * r1 ** 2 * r2 ** 2 * sin(theta) ** 2
           * (sinphi1 * cos(phi - phi2) - sinphi2 * cos(phi - phi1)) ** 2
           / c ** 2)

        A = (r1 ** 2 * sinphi1 ** 2
             + r2 ** 2 * sinphi2 ** 2
             - r1 * r2 * sinphi1 * sinphi2)
        B = (2 * r1 ** 2 * sinphi1 * cosphi1
             + 2 * r2 ** 2 * sinphi2 * cosphi2
             - r1 * r2 * sinphi2 * cosphi1
             - r1 * r2 * sinphi1 * cosphi2)
        C = (r1 ** 2 * cosphi1 ** 2
             + r2 ** 2 * cosphi2 ** 2
             - r1 * r2 * cosphi1 * cosphi2)

        return 2 * (A * tanphi ** 2 + B * tanphi + C) / den

    def dphi_dt0(self, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta)
               * (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2))
               / c)
        num = (r2 * cosphi2 - r1 * cosphi1
               + tanphi * (r2 * sinphi2 - r1 * sinphi1))

        return num / den

    def dphi_dt1(self, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta)
               * (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2))
               / c)
        num = -r2 * (sinphi2 * tanphi + cosphi2)

        return num / den

    def dphi_dt2(self, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta)
               * (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2))
               / c)
        num = r1 * (sinphi1 * tanphi + cosphi1)

        return num / den

def do_full_reconstruction(data, tablename):
    """Do a reconstruction of all simulated data and store results"""

    try:
        data.createGroup('/', 'reconstructions', "Angle reconstructions")
    except tables.NodeError:
        pass

    try:
        data.removeNode('/reconstructions', tablename)
    except tables.NoSuchNodeError:
        pass

    table = data.createTable('/reconstructions', tablename,
                             ReconstructedEvent, "Reconstruction data")

    rec = DirectionReconstruction(data, table, min_n134=1, N=100)

    rec.reconstruct_angles(tablename='zenith_0', THETA=0)
    rec.reconstruct_angles(tablename='zenith_5', THETA=deg2rad(5))
    rec.reconstruct_angles(tablename='zenith_22_5', THETA=pi / 8)
    rec.reconstruct_angles(tablename='zenith_35', THETA=deg2rad(35))

    # SPECIALS
    # Station sizes
    rec.reconstruct_angles(tablename='zenith_22_5_size5', THETA=pi / 8)
    rec.reconstruct_angles(tablename='zenith_22_5_size20', THETA=pi / 8)

    # SPECIALS
    # Binnings
    kwargs = dict(tablename='zenith_22_5',
                  THETA=pi / 8)
    kwargs['randomize_binning'] = False
    rec.reconstruct_angles(binning=1, **kwargs)
    rec.reconstruct_angles(binning=2.5, **kwargs)
    rec.reconstruct_angles(binning=5, **kwargs)
    kwargs['randomize_binning'] = True
    rec.reconstruct_angles(binning=1, **kwargs)
    rec.reconstruct_angles(binning=2.5, **kwargs)
    rec.reconstruct_angles(binning=5, **kwargs)

def do_reconstruction_plots(data, tablename):
    """Make plots based upon earlier reconstructions"""

    table = data.getNode('/reconstructions', tablename)

    plot_uncertainty_mip(table)
    plot_uncertainty_zenith(table)
    plot_uncertainty_size(table)
    plot_uncertainty_binsize(table)

    plot_phi_reconstruction_results_for_MIP(table, 1)
    plot_phi_reconstruction_results_for_MIP(table, 2)
    boxplot_theta_reconstruction_results_for_MIP(table, 1)
    boxplot_theta_reconstruction_results_for_MIP(table, 2)
    boxplot_phi_reconstruction_results_for_MIP(table, 1)
    boxplot_phi_reconstruction_results_for_MIP(table, 2)
    boxplot_arrival_times(table, 1)
    boxplot_arrival_times(table, 2)
    boxplot_core_distances_for_mips(table)
    plot_detection_efficiency_vs_R_for_angles(1)
    plot_detection_efficiency_vs_R_for_angles(2)
    plot_reconstruction_efficiency_vs_R_for_angles(1)
    plot_reconstruction_efficiency_vs_R_for_angles(2)
    plot_reconstruction_efficiency_vs_R_for_mips()

def plot_uncertainty_mip(table):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    figure()
    rcParams['text.usetex'] = False
    x, y, y2 = [], [], []
    for D in range(1, 6):
        x.append(D)
        events = table.readWhere(
            '(min_n134==%d) & (reference_theta==%.40f) & (size==10) & (bin==0) & '
            '(0 < r) & (r <= 100)' % (D, float32(pi / 8)))
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
    # Uncertainty estimate
    x = linspace(1, 5, 50)
    phis = linspace(-pi, pi, 50)
    phi_errsq = mean(rel_phi_errorsq(pi / 8, phis, phi1, phi2))
    theta_errsq = mean(rel_theta_errorsq(pi / 8, phis, phi1, phi2))
    y = TIMING_ERROR * std_t(x) * sqrt(phi_errsq)
    y2 = TIMING_ERROR * std_t(x) * sqrt(theta_errsq)
    plot(x, rad2deg(y), label="Estimate Phi")
    plot(x, rad2deg(y2), label="Estimate Theta")
    # Labels etc.
    xlabel("Minimum number of particles")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$\theta = 22.5^\circ$")
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    utils.saveplot()
    print

def plot_uncertainty_zenith(table):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    figure()
    rcParams['text.usetex'] = False
    x, y, y2 = [], [], []
    global MYT
    MYT = []
    for THETA in [0, deg2rad(5), pi / 8, deg2rad(35)]:
        x.append(THETA)
        events = table.readWhere(
            '(min_n134==2) & (reference_theta==%.40f) & (size==10) & (bin==0)' %
            float32(THETA))
        MYT.append((events['t1'], events['t3'], events['t4']))
        print rad2deg(THETA), len(events),
        errors = events['reference_theta'] - events['reconstructed_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['reference_phi'] - events['reconstructed_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(rad2deg(x), rad2deg(y), '^', label="Theta")
    # Azimuthal angle undefined for zenith = 0
    plot(rad2deg(x[1:]), rad2deg(y2[1:]), 'v', label="Phi")
    print
    print "zenith: theta, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    # Uncertainty estimate
    x = linspace(0, deg2rad(35), 50)
    phis = linspace(-pi, pi, 50)
    y, y2, y3 = [], [], []
    for t in x:
        y.append(mean(rel_phi_errorsq(t, phis, phi1, phi2)))
        y3.append(mean(rel_phi_errorsq(t, phis, phi1, phi2)) * sin(t) ** 2)
        y2.append(mean(rel_theta_errorsq(t, phis, phi1, phi2)))
    y = TIMING_ERROR * sqrt(array(y))
    y3 = TIMING_ERROR * sqrt(array(y3))
    y2 = TIMING_ERROR * sqrt(array(y2))
    plot(rad2deg(x), rad2deg(y), label="Estimate Phi")
    plot(rad2deg(x), rad2deg(y3), label="Estimate Phi * sin(Theta)")
    plot(rad2deg(x), rad2deg(y2), label="Estimate Theta")
    # Labels etc.
    xlabel("Shower zenith angle (degrees)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$N_{MIP} = 2$")
    ylim(0, 100)
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    utils.saveplot()
    print

def plot_uncertainty_size(table):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    figure()
    rcParams['text.usetex'] = False
    x, y, y2 = [], [], []
    for size in [5, 10, 20]:
        x.append(size)
        events = table.readWhere(
            '(min_n134==2) & (reference_theta==%.40f) & (size==%d) & (bin==0)' %
            (float32(pi / 8), size))
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
    print "stationsize: size, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    # Uncertainty estimate
    x = linspace(5, 20, 50)
    phis = linspace(-pi, pi, 50)
    y, y2 = [], []
    for s in x:
        y.append(mean(rel_phi_errorsq(pi / 8, phis, phi1, phi2, r1=s, r2=s)))
        y2.append(mean(rel_theta_errorsq(pi / 8, phis, phi1, phi2, r1=s, r2=s)))
    y = TIMING_ERROR * sqrt(array(y))
    y2 = TIMING_ERROR * sqrt(array(y2))
    plot(x, rad2deg(y), label="Estimate Phi")
    plot(x, rad2deg(y2), label="Estimate Theta")
    # Labels etc.
    xlabel("Station size (m)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$\theta = 22.5^\circ, N_{MIP} = 2$")
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    utils.saveplot()
    print

def plot_uncertainty_binsize(table):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    figure()
    rcParams['text.usetex'] = False
    x, y, y2 = [], [], []
    for bin_size in [0, 1, 2.5, 5]:
        if bin_size == 0:
            is_randomized = False
        else:
            is_randomized = True
        x.append(bin_size)
        events = table.readWhere(
            '(min_n134==2) & (reference_theta==%.40f) & (size==10) & (bin==%.40f) & '
            '(bin_r==%s)' %
            (float32(pi / 8), bin_size, is_randomized))
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
    print "binsize: size, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    # Uncertainty estimate
    x = linspace(0, 5, 50)
    phis = linspace(-pi, pi, 50)
    y, y2 = [], []
    phi_errorsq = mean(rel_phi_errorsq(pi / 8, phis, phi1, phi2))
    theta_errorsq = mean(rel_theta_errorsq(pi / 8, phis, phi1, phi2))
    for t in x:
        y.append(sqrt((TIMING_ERROR ** 2 + t ** 2 / 12) * phi_errorsq))
        y2.append(sqrt((TIMING_ERROR ** 2 + t ** 2 / 12) * theta_errorsq))
    y = array(y)
    y2 = array(y2)
    plot(x, rad2deg(y), label="Estimate Phi")
    plot(x, rad2deg(y2), label="Estimate Theta")
    # Labels etc.
    xlabel("Bin size (ns)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$\theta = 22.5^\circ, N_{MIP} = 2$")
    legend(loc='best', numpoints=1)
    ylim(ymin=0)
    if USE_TEX:
        rcParams['text.usetex'] = True
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

def plot_phi_reconstruction_results_for_MIP(table, N):
    min_theta = 1
    query = '(min_n134>=%d) & (size==10) & (bin==0) & (reference_theta > %.5f)' % \
            (N, deg2rad(min_theta))

    events = table.readWhere(query)
    sim_phi = events['reference_phi']
    r_phi = events['reconstructed_phi']

    figure()
    plot_2d_histogram(rad2deg(sim_phi), rad2deg(r_phi), 180)
    xlabel(r"$\phi_{simulated}$")
    ylabel(r"$\phi_{reconstructed}$")
    title(r"$N_{MIP} \geq %d, \quad \theta > 0^\circ$" % N)

    utils.saveplot(N)

def boxplot_theta_reconstruction_results_for_MIP(table, N):
    figure()
    query = '(min_n134>=%d) & (size==10) & (bin==0)' % N

    angles = [0, 5, 22.5, 35]
    r_dtheta = []
    for angle in angles:
        low, high = deg2rad(angle - 1), deg2rad(angle + 1)
        sel_query = query + '& (low < reference_theta) & (reference_theta < high)'
        sel = table.readWhere(sel_query)
        r_dtheta.append(rad2deg(sel[:]['reconstructed_theta'] - sel[:]['reference_theta']))

    boxplot(r_dtheta, sym='')
    xticks(range(1, len(angles) + 1), [str(u) for u in angles])

    xlabel(r"$\theta_{simulated}$")
    ylabel(r"$\theta_{reconstructed} - \theta_{simulated}$")
    title(r"$N_{MIP} \geq %d$" % N)

    axhline(0)
    ylim(-20, 20)

    utils.saveplot(N)

def boxplot_phi_reconstruction_results_for_MIP(table, N):
    figure()
    min_theta = 1
    query = '(min_n134>=%d) & (size==10) & (bin==0) & (reference_theta > %.5f)' % \
            (N, deg2rad(min_theta))

    bin_edges = linspace(-180, 180, 18)
    x, r_dphi = [], []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        rad_low = deg2rad(low)
        rad_high = deg2rad(high)
        sel_query = query + '& (rad_low < reference_phi) & (reference_phi < rad_high)'
        sel = table.readWhere(sel_query)
        dphi = sel[:]['reconstructed_phi'] - sel[:]['reference_phi']
        dphi = (dphi + pi) % (2 * pi) - pi
        r_dphi.append(rad2deg(dphi))
        x.append((low + high) / 2)

    boxplot(r_dphi, positions=x, widths=.7 * (high - low), sym='')

    xlabel(r"$\phi_{simulated}$")
    ylabel(r"$\phi_{reconstructed} - \phi_{simulated}$")
    title(r"$N_{MIP} \geq %d, \quad \theta > 0^\circ$" % N)

    xticks(linspace(-180, 180, 9))
    axhline(0)

    utils.saveplot(N)

def boxplot_arrival_times(table, N):
    figure()
    query = '(min_n134>=%d) & (size==10) & (bin==0) & (reference_theta == 0)' % N

    bin_edges = linspace(0, 100, 10)
    x, arrival_times = [], []
    for low, high in zip(bin_edges[:-1], bin_edges[1:]):
        sel_query = query + '& (low <= r) & (r < high)'
        sel = table.readWhere(sel_query)
        t2 = sel[:]['t2']
        arrival_times.append(t2.compress(t2 > -999))
        x.append((low + high) / 2)

    boxplot(arrival_times, positions=x, widths=.7 * (high - low), sym='')

    xlabel("Core distance [m]")
    ylabel("Arrival time [ns]")
    title(r"$N_{MIP} \geq %d, \quad \theta = 0^\circ$" % N)

    xticks(arange(0, 100.5, 10))

    utils.saveplot(N)

def boxplot_core_distances_for_mips(table):
    query = ('(size==10) & (bin==0) & (reference_theta == %.40f)' %
             float32(deg2rad(22.5)))

    r_list = []
    for D in range(1, 5):
        sel_query = query + '& (min_n134 == %d)' % D
        sel = table.readWhere(sel_query)
        r = sel[:]['r']
        r_list.append(r)

    figure()
    boxplot(r_list)
    xticks(range(1, 5))
    xlabel("Minimum number of particles")
    ylabel("Core distance [m]")
    title(r"$\theta = 22.5^\circ$")

    utils.saveplot()

def plot_detection_efficiency_vs_R_for_angles(N):
    figure()

    bin_edges = linspace(0, 100, 20)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2.

    for angle in [0, 5, 22.5, 35]:
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
    title(r"$N_{MIP} \geq %d$" % N)
    legend()

    utils.saveplot(N)

def plot_reconstruction_efficiency_vs_R_for_angles(N):
    figure()

    bin_edges = linspace(0, 100, 10)
    x = (bin_edges[:-1] + bin_edges[1:]) / 2.

    for angle in [0, 5, 22.5, 35]:
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
                shower_results.append(len(sel))
            query = "(size == 10) & (bin == 0) & (reference_theta == %.40f) & (low <= r) & (r < high) & (min_n134 >= N)" % float32(deg2rad(angle))
            ssel = data.root.reconstructions.full.readWhere(query)
            efficiencies.append(len(ssel) / sum(shower_results))

        plot(x, efficiencies, label=r'$\theta = %s^\circ$' % angle)

    xlabel("Core distance [m]")
    ylabel("Reconstruction efficiency")
    title(r"$N_{MIP} \geq %d$" % N)
    legend()

    utils.saveplot(N)

def plot_reconstruction_efficiency_vs_R_for_mips():
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
            query = "(size == 10) & (bin == 0) & (reference_theta == %.40f) & (low <= r) & (r < high) & (min_n134 == N)" % float32(deg2rad(22.5))
            ssel = data.root.reconstructions.full.readWhere(query)
            print sum(shower_results), len(ssel), len(ssel) / sum(shower_results)
            efficiencies.append(len(ssel) / sum(shower_results))

        plot(x, efficiencies, label=r'$N_{MIP} = %d' % N)

    xlabel("Core distance [m]")
    ylabel("Reconstruction efficiency")
    title(r"$\theta = 22.5^\circ$")
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
        data = tables.openFile('master.h5', 'a')

#    if '/reconstructions/full' not in data:
#        print "Reconstructing shower direction..."
#        do_full_reconstruction(data, 'full')
#    else:
#        print "Skipping reconstruction!"

    print "Reconstructing shower direction..."
    do_full_reconstruction(data, 'fulltest')

    #utils.set_prefix("DIR-")
    #do_reconstruction_plots(data, 'full')

#    utils.set_prefix("WIP-")
#    do_jos_plots(data)
