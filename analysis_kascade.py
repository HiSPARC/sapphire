import tables
from itertools import combinations
import re
import csv
import zlib

from pylab import *
from scipy.optimize import curve_fit

from scipy import integrate
from scipy.special import erf

from itertools import combinations

from hisparc.analysis import kascade_coincidences
from hisparc.containers import Coincidence

DATAFILE = 'generator.h5'

USE_TEX = False


ADC_THRESHOLD = 20
ADC_TIME_PER_SAMPLE = 2.5e-9
ADC_MIP = 400.

D_Z = 1 # Delta zenith (used in data selection)

DETECTORS = [(65., 15.05, 'UD'), (65., 20.82, 'UD'), (70., 23.71, 'LR'),
             (60., 23.71, 'LR')]

TIMING_ERROR = 4
#TIMING_ERROR = 7

# For matplotlib plots
if USE_TEX:
    rcParams['font.serif'] = 'Computer Modern'
    rcParams['font.sans-serif'] = 'Computer Modern'
    rcParams['font.family'] = 'sans-serif'
    rcParams['figure.figsize'] = [5 * x for x in (1, 2./3)]
    rcParams['figure.subplot.left'] = 0.125
    rcParams['figure.subplot.bottom'] = 0.125
    rcParams['font.size'] = 11
    rcParams['legend.fontsize'] = 'small'


def process_traces(events, traces_table, limit=None):
    """Process traces to yield pulse timing information"""

    result = []
    result2 = []
    for event in events[:limit]:
        trace = get_traces(traces_table, event)
        timings, timings2 = zip(*[reconstruct_time_from_trace(x) for x in
                                  trace])
        result.append(timings)
        result2.append(timings2)
    return 1e9 * array(result), 1e9 * array(result2)

def get_traces(traces_table, event):
    """Retrieve traces from table and reconstruct them"""

    if type(event) != list:
        idxs = event['traces']
    else:
        idxs = event

    traces = []
    for idx in idxs:
        trace = zlib.decompress(traces_table[idx]).split(',')
        if trace[-1] == '':
            del trace[-1]
        trace = array([int(x) for x in trace])
        traces.append(trace)
    return traces

def reconstruct_time_from_trace(trace):
    """Reconstruct time of measurement from a trace"""

    t = trace[:100]
    baseline = mean(t)
    stdev = std(t)

    trace = trace - baseline
    threshold = ADC_THRESHOLD

    value = nan
    for i, t in enumerate(trace):
        if t >= threshold:
            value = i
            break

    # Better value, interpolation
    if not isnan(value):
        x0, x1 = i - 1, i
        y0, y1 = trace[x0], trace[x1]
        v2 = 1. * (threshold - y0) / (y1 - y0) + x0
    else:
        v2 = nan

    return value * ADC_TIME_PER_SAMPLE, v2 * ADC_TIME_PER_SAMPLE

class ReconstructedEvent(tables.IsDescription):
    r = tables.Float32Col()
    t1 = tables.Float32Col()
    t2 = tables.Float32Col()
    t3 = tables.Float32Col()
    t4 = tables.Float32Col()
    n1 = tables.UInt16Col()
    n2 = tables.UInt16Col()
    n3 = tables.UInt16Col()
    n4 = tables.UInt16Col()
    k_theta = tables.Float32Col()
    k_phi = tables.Float32Col()
    k_energy = tables.Float32Col()
    h_theta = tables.Float32Col()
    h_theta1 = tables.Float32Col()
    h_theta2 = tables.Float32Col()
    h_phi = tables.Float32Col()
    D = tables.UInt16Col()

def reconstruct_angle(event, R=10):
    """Reconstruct angles from a single event"""

    dt1 = event['t1'] - event['t3']
    dt2 = event['t1'] - event['t4']

    return reconstruct_angle_dt(dt1, dt2, R)

def reconstruct_angle_dt(dt1, dt2, R=10):
    """Reconstruct angle given time differences"""

    c = 3.00e+8

    r1 = R
    r2 = R 

    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                  (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
    theta1 = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
    theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))

    e1 = sqrt(rel_theta1_errorsq(theta1, phi, phi1, phi2, R, R))
    e2 = sqrt(rel_theta2_errorsq(theta2, phi, phi1, phi2, R, R))

    theta_wgt = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

    return theta_wgt, phi, theta1, theta2

def calc_phi(s1, s2):
    """Calculate angle between detectors (phi1, phi2)"""

    x1, y1 = DETECTORS[s1 - 1][:2]
    x2, y2 = DETECTORS[s2 - 1][:2]

    return arctan2((y2 - y1), (x2 - x1))

def rel_phi_errorsq(theta, phi, phi1, phi2, r1=10, r2=10):
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

def dphi_dt0(theta, phi, phi1, phi2, r1=10, r2=10):
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

def dphi_dt1(theta, phi, phi1, phi2, r1=10, r2=10):
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

def dphi_dt2(theta, phi, phi1, phi2, r1=10, r2=10):
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

def rel_theta_errorsq(theta, phi, phi1, phi2, r1=10, r2=10):
    e1 = rel_theta1_errorsq(theta, phi, phi1, phi2, r1, r2)
    e2 = rel_theta2_errorsq(theta, phi, phi1, phi2, r1, r2)

    #return minimum(e1, e2)
    return e1

def rel_theta1_errorsq(theta, phi, phi1, phi2, r1=10, r2=10):
    # speed of light in m / ns
    c = .3

    sintheta = sin(theta)
    sinphiphi1 = sin(phi - phi1)

    den = (1 - sintheta ** 2) * r1 ** 2 * cos(phi - phi1) ** 2

    A = (r1 ** 2 * sinphiphi1 ** 2
         * rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
    B = (r1 * c * sinphiphi1
         * (dphi_dt0(theta, phi, phi1, phi2, r1, r2)
            - dphi_dt1(theta, phi, phi1, phi2, r1, r2)))
    C = 2 * c ** 2

    errsq = (A * sintheta ** 2 - 2 * B * sintheta + C) / den

    return where(isnan(errsq), inf, errsq)

def rel_theta2_errorsq(theta, phi, phi1, phi2, r1=10, r2=10):
    # speed of light in m / ns
    c = .3

    sintheta = sin(theta)
    sinphiphi2 = sin(phi - phi2)

    den = (1 - sintheta ** 2) * r2 ** 2 * cos(phi - phi2) ** 2

    A = (r2 ** 2 * sinphiphi2 ** 2
         * rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
    B = (r2 * c * sinphiphi2
         * (dphi_dt0(theta, phi, phi1, phi2, r1, r2)
            - dphi_dt2(theta, phi, phi1, phi2, r1, r2)))
    C = 2 * c ** 2

    errsq = (A * sintheta ** 2 - 2 * B * sintheta + C) / den

    return where(isnan(errsq), inf, errsq)

def reconstruct_angles(data, dstname, events, timing_data, shifts=None,
                       N=None):
    """Reconstruct angles"""

    try:
        data.createGroup('/', 'reconstructions', "Angle reconstructions")
    except tables.NodeError:
        pass

    try:
        data.removeNode('/reconstructions', dstname)
    except tables.NoSuchNodeError:
        pass

    dest = data.createTable('/reconstructions', dstname,
                            ReconstructedEvent, "Reconstruction data")

    R = 10

    NT, NS = 0, 0

    dst_row = dest.row
    for rawevent, timing in zip(events[:N], timing_data):
        NT += 1
        n1, n2, n3, n4 = rawevent['pulseheights'] / ADC_MIP
        if shifts:
            timing += shifts

        event = dict(n1=n1, n2=n2, n3=n3, n4=n4, t1=timing[0],
                     t2=timing[1], t3=timing[2], t4=timing[3])

        theta, phi, theta1, theta2 = reconstruct_angle(event, R)

        if not isnan(theta) and not isnan(phi):
            NS += 1
            xc, yc = rawevent['k_core_pos'] - DETECTORS[1][:2]
            dst_row['r'] = sqrt(xc ** 2 + yc ** 2)
            dst_row['t1'] = event['t1']
            dst_row['t2'] = event['t2']
            dst_row['t3'] = event['t3']
            dst_row['t4'] = event['t4']
            dst_row['n1'] = event['n1']
            dst_row['n2'] = event['n2']
            dst_row['n3'] = event['n3']
            dst_row['n4'] = event['n4']
            dst_row['k_theta'] = rawevent['k_zenith']
            dst_row['k_phi'] = -(rawevent['k_azimuth'] + deg2rad(75)) % \
                               (2 * pi) - pi
            dst_row['k_energy'] = rawevent['k_energy']
            dst_row['h_theta'] = theta
            dst_row['h_theta1'] = theta1
            dst_row['h_theta2'] = theta2
            dst_row['h_phi'] = phi
            dst_row['D'] = round(min(event['n1'], event['n3'], event['n4']))
            dst_row.append()
    dest.flush()

    print NT, NS

def do_reconstruction_plots(data, tablename, table2name, sim_data,
                            sim_tablename):
    """Make plots based upon earlier reconstructions"""

    table = data.getNode('/reconstructions', tablename)
    table2 = data.getNode('/reconstructions', table2name)
    sim_table = sim_data.getNode('/reconstructions', sim_tablename)

    plot_2d_results_phi(data)
    plot_2d_results_theta(data)

    plot_uncertainty_mip(table, sim_table)
    plot_uncertainty_zenith(table, sim_table)
    plot_uncertainty_zenith2(table, table2)
    plot_uncertainty_energy(table)
    plot_mip_core_dists_mean(table, sim_table)
    plot_zenith_core_dists_mean(table, sim_table)
    plot_uncertainty_core_dist_phi_theta(table, sim_table)

def plot_uncertainty_mip(table, sim_table):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    THETA = pi / 8

    figure()
    rcParams['text.usetex'] = False
    x, y, y2, sy, sy2 = [], [], [], [], []
    for D in range(1, 6):
        x.append(D)
        ### KASCADE data
        events = table.readWhere(
            '(D==%d) & (%f <= k_theta) & (k_theta < %f)'
             % (D, THETA - deg2rad(D_Z), THETA + deg2rad(D_Z)))
        print len(events),
        errors = events['k_theta'] - events['h_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['k_phi'] - events['h_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))

        ### simulation data
        events = sim_table.readWhere(
            '(D==%d) & (sim_theta==%.40f) & (size==10) & (bin==0) & '
            '(0 < r) & (r <= 100)' % (D, float32(THETA)))
        print len(events),
        errors = events['sim_theta'] - events['r_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['sim_phi'] - events['r_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        sy.append(std(errors))
        sy2.append(std(errors2))
    print
    print "mip: D, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    # Uncertainty estimate
    cx = linspace(1, 5, 50)
    phis = linspace(-pi, pi, 50)
    phi_errsq = mean(rel_phi_errorsq(pi / 8, phis, phi1, phi2))
    theta_errsq = mean(rel_theta_errorsq(pi / 8, phis, phi1, phi2))
    cy = TIMING_ERROR * std_t(cx) * sqrt(phi_errsq)
    cy2 = TIMING_ERROR * std_t(cx) * sqrt(theta_errsq)

    plot(x, rad2deg(y), '^', label="Theta")
    plot(x, rad2deg(sy), '^', label="Theta (simulation)")
    plot(cx, rad2deg(cy2), label="Theta (calculation)")
    plot(x, rad2deg(y2), 'v', label="Phi")
    plot(x, rad2deg(sy2), 'v', label="Phi (simulation)")
    plot(cx, rad2deg(cy), label="Phi (calculation)")
    # Labels etc.
    xlabel("$N_{MIP}$")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$\theta = %.1f \pm %d^\circ$" % (rad2deg(THETA), D_Z))
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-MIP.pdf')
    print

def plot_uncertainty_zenith(table, sim_table):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    figure()
    rcParams['text.usetex'] = False
    x, y, y2, sy, sy2 = [], [], [], [], []
    for THETA in [0, deg2rad(5), pi / 8, deg2rad(35)]:
        x.append(THETA)
        ### KASCADE data
        events = table.readWhere(
            '(D==2) & (%f <= k_theta) & (k_theta < %f)'
             % (THETA - deg2rad(D_Z), THETA + deg2rad(D_Z)))
        print rad2deg(THETA), len(events),
        errors = events['k_theta'] - events['h_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['k_phi'] - events['h_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))

        ### simulation data
        events = sim_table.readWhere(
            '(D==2) & (sim_theta==%.40f) & (size==10) & (bin==0)' %
            float32(THETA))
        print rad2deg(THETA), len(events),
        errors = events['sim_theta'] - events['r_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['sim_phi'] - events['r_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        sy.append(std(errors))
        sy2.append(std(errors2))
    print
    print "zenith: theta, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    # Uncertainty estimate
    cx = linspace(0, deg2rad(35), 50)
    phis = linspace(-pi, pi, 50)
    cy, cy2, cy3 = [], [], []
    for t in cx:
        cy.append(mean(rel_phi_errorsq(t, phis, phi1, phi2)))
        cy3.append(mean(rel_phi_errorsq(t, phis, phi1, phi2)) * sin(t) ** 2)
        cy2.append(mean(rel_theta_errorsq(t, phis, phi1, phi2)))
    cy = TIMING_ERROR * sqrt(array(cy))
    cy3 = TIMING_ERROR * sqrt(array(cy3))
    cy2 = TIMING_ERROR * sqrt(array(cy2))


    plot(rad2deg(x), rad2deg(y), '^', label="Theta")
    plot(rad2deg(x), rad2deg(sy), '^', label="Theta (simulation)")
    plot(rad2deg(cx), rad2deg(cy2), label="Theta (calculation)")
    # Azimuthal angle undefined for zenith = 0
    plot(rad2deg(x[1:]), rad2deg(y2[1:]), 'v', label="Phi")
    plot(rad2deg(x[1:]), rad2deg(sy2[1:]), 'v', label="Phi (simulation)")
    plot(rad2deg(cx), rad2deg(cy), label="Phi (calculation)")
    plot(rad2deg(cx), rad2deg(cy3), label="Phi (calculation) * sin(Theta)")
    # Labels etc.
    xlabel("Zenith angle (deg $\pm %d$)" % D_Z)
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$N_{MIP} = 2$")
    ylim(0, 100)
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-zenith.pdf')
    print

def plot_uncertainty_zenith2(table, table2):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    figure()
    rcParams['text.usetex'] = False
    x, y, y2, ty, ty2 = [], [], [], [], []
    for THETA in [deg2rad(5), pi / 8, deg2rad(35)]:
        x.append(THETA)
        ### KASCADE data
        events = table.readWhere(
            '(D==2) & (%f <= k_theta) & (k_theta < %f)'
             % (THETA - deg2rad(D_Z), THETA + deg2rad(D_Z)))
        print rad2deg(THETA), len(events),
        errors = events['k_theta'] - events['h_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['k_phi'] - events['h_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2) * sin(THETA))

        ### other KASCADE data
        events = table2.readWhere(
            '(D==2) & (%f <= k_theta) & (k_theta < %f)'
             % (THETA - deg2rad(D_Z), THETA + deg2rad(D_Z)))
        print rad2deg(THETA), len(events),
        errors = events['k_theta'] - events['h_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['k_phi'] - events['h_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        ty.append(std(errors))
        ty2.append(std(errors2) * sin(THETA))
    print
    print "zenith: theta, theta_std, phi_std"
    for u, v, w in zip(x, y, y2):
        print u, v, w
    print
    
    errorbar(rad2deg(x), rad2deg(y), xerr=D_Z, fmt='^', label="Theta (FSOT)")
    errorbar(rad2deg(x), rad2deg(ty), xerr=D_Z, fmt='^', label="Theta (LINT)")
    # Azimuthal angle undefined for zenith = 0
    errorbar(rad2deg(x[1:]), rad2deg(y2[1:]), xerr=D_Z, fmt='v', label="Phi (FSOT)")
    errorbar(rad2deg(x[1:]), rad2deg(ty2[1:]), xerr=D_Z, fmt='v', label="Phi (LINT)")
    # Labels etc.
    xlabel("Zenith angle (deg)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$N_{MIP} = 2$")
    legend(loc='lower right', numpoints=1)
    ylim(0, 14)
    xlim(0, 40)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-zenith2.pdf')
    print

def plot_uncertainty_phi(table):
    # constants for uncertainty estimation
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    figure()
    rcParams['text.usetex'] = False
    # Uncertainty estimate
    x = linspace(0, deg2rad(360), 50)
    y, y2 = [], []
    for p in x:
        y.append(rel_phi_errorsq(pi / 8, p, phi1, phi2))
        y2.append(rel_theta_errorsq(pi / 8, p, phi1, phi2))
    y = TIMING_ERROR * sqrt(array(y))
    y2 = TIMING_ERROR * sqrt(array(y2))
    plot(rad2deg(x), rad2deg(y), label="Estimate Phi")
    plot(rad2deg(x), rad2deg(y2), label="Estimate Theta")
    # Labels etc.
    xlabel("Shower azimuth angle (degrees)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$\theta = 22.5^\circ, N_{MIP} = 2$")
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-phi.pdf')
    print

def plot_uncertainty_energy(table):
    THETA = pi / 8
    N = 2
    D2_Z = 5 * D_Z
    D_lE = .2

    figure()
    rcParams['text.usetex'] = False
    x, y, y2 = [], [], []
    for lE in range(14, 18):
        x.append(lE)
        ### KASCADE data
        events = table.readWhere(
            '(D==%d) & (%f <= k_theta) & (k_theta < %f) & '
            '(%f <= k_energy) & (k_energy < %f)'
             % (N, THETA - deg2rad(D2_Z), THETA + deg2rad(D2_Z),
                10 ** (lE - D_lE), 10 ** (lE + D_lE)))
        print len(events),
        errors = events['k_theta'] - events['h_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['k_phi'] - events['h_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))

    plot(x, rad2deg(y), '^', label="Theta")
    plot(x, rad2deg(y2), 'v', label="Phi")
    # Labels etc.
    xlabel("log Energy (eV) $\pm %.1f$" % D_lE)
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$\theta = %.1f \pm %d^\circ$" % (rad2deg(THETA), D2_Z))
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-energy.pdf')
    print

# Time of first hit pamflet functions
Q = lambda t, n: ((.5 * (1 - erf(t / sqrt(2)))) ** (n - 1)
                  * exp(-.5 * t ** 2) / sqrt(2 * pi))

expv_t = vectorize(lambda n: integrate.quad(lambda t: t * Q(t, n)
                                                      / n ** -1,
                                            -inf, +inf))
expv_tv = lambda n: expv_t(n)[0]
expv_tsq = vectorize(lambda n: integrate.quad(lambda t: t ** 2 * Q(t, n)
                                                        / n ** -1,
                                              -inf, +inf))
expv_tsqv = lambda n: expv_tsq(n)[0]

std_t = lambda n: sqrt(expv_tsqv(n) - expv_tv(n) ** 2)

def plot_mip_core_dists_mean(table, sim_table):
    figure()
    rcParams['text.usetex'] = False
    THETA = pi / 8
    x, y, yerr, sy, syerr = [], [], [], [], []
    for N in range(1, 6):
        events = table.readWhere(
            '(D==%d) & (%f <= k_theta) & (k_theta < %f) & '
            '(%f <= k_energy) & (k_energy < %f)'
             % (N, THETA - deg2rad(D_Z), THETA + deg2rad(D_Z),
                8e14, 2e15))
        x.append(N)
        y.append(mean(events['r']))
        yerr.append(std(events['r']))

        events = sim_table.readWhere(
            '(D==%d) & (sim_theta==%.40f) & (size==10) & (bin==0)' %
            (N, float32(THETA)))
        sy.append(mean(events['r']))
        syerr.append(std(events['r']))
    plot(x, y, '^-', label=r"$(.8\,\mathrm{PeV} \leq E < 2\,\mathrm{PeV})$")
    plot(x, sy, '^-', label=r"(simulation, $E = 1\,\mathrm{PeV}$)")
    print "zenith: theta, r_mean, r_std"
    for u, v, w in zip(x, y, yerr):
        print u, v, w
    print
    # Labels etc.
    xlabel("$N_{MIP}$")
    ylabel("Core distance (m)")
    title(r"$\theta = %.1f \pm %d^\circ$" % (rad2deg(THETA), D_Z))
    xlim(.5, 5.5)
    ylim(ymin=0)
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-mip_core_dists_mean.pdf')
    print

def plot_zenith_core_dists_mean(table, sim_table):
    figure()
    rcParams['text.usetex'] = False
    x, y, yerr, sy, syerr = [], [], [], [], []
    N = 2
    for THETA in [0, deg2rad(5), pi / 8, deg2rad(35)]:
        events = table.readWhere(
            '(D==%d) & (%f <= k_theta) & (k_theta < %f) & '
            '(%f <= k_energy) & (k_energy < %f)'
             % (N, THETA - deg2rad(D_Z), THETA + deg2rad(D_Z),
                8e14, 2e15))
        x.append(THETA)
        y.append(mean(events['r']))
        yerr.append(std(events['r']))

        events = sim_table.readWhere(
            '(D==%d) & (sim_theta==%.40f) & (size==10) & (bin==0)' %
            (N, float32(THETA)))
        sy.append(mean(events['r']))
        syerr.append(std(events['r']))
    plot(rad2deg(x), y, '^-',
         label=r"$(.8\,\mathrm{PeV} \leq E < 2\,\mathrm{PeV})$")
    plot(rad2deg(x), sy, '^-', label=r"(simulation, $E = 1\,\mathrm{PeV}$)")
    # Labels etc.
    xlabel("Zenith angle (deg $\pm %d$)" % D_Z)
    ylabel("Core distance (m)")
    title(r"$N_{MIP} = %d$" % N)
    legend(numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-zenith_core_dists_mean.pdf')

def plot_uncertainty_core_dist_phi_theta(table, sim_table):
    THETA = pi / 8
    bins = linspace(0, 80, 6)
    N = 2
    D2_Z = 5
    events = table.readWhere(
        '(D==2) & (%f <= k_theta) & (k_theta < %f) & '
        '(%f <= k_energy) & (k_energy < %f)'
         % (THETA - deg2rad(D2_Z), THETA + deg2rad(D2_Z),
            8e14, 2e15))
    sim_events = sim_table.readWhere(
            '(D==2) & (sim_theta==%.40f) & (size==10) & (bin==0)' %
            (float32(THETA)))
    x, y, y2, l, sx, sy, sy2 = [], [], [], [], [], [], []
    for r0, r1 in zip(bins[:-1], bins[1:]):
        e = events.compress((r0 <= events['r']) & (events['r'] < r1))
        if len(e) > 10:
            errors = e['k_phi'] - e['h_phi']
            errors2 = e['k_theta'] - e['h_theta']
            # Make sure -pi < errors < pi
            errors = (errors + pi) % (2 * pi) - pi
            errors2 = (errors2 + pi) % (2 * pi) - pi
            x.append(mean([r0, r1]))
            y.append(std(errors))
            y2.append(std(errors2))
            l.append(len(e))

        e = sim_events.compress((r0 <= sim_events['r']) & 
                                (sim_events['r'] < r1))
        if len(e) > 10:
            errors = e['sim_phi'] - e['r_phi']
            errors2 = e['sim_theta'] - e['r_theta']
            # Make sure -pi < errors < pi
            errors = (errors + pi) % (2 * pi) - pi
            errors2 = (errors2 + pi) % (2 * pi) - pi
            sx.append(mean([r0, r1]))
            sy.append(std(errors))
            sy2.append(std(errors2))

    figure()
    rcParams['text.usetex'] = False
    plot(x, rad2deg(y2), '^-',
         label="Theta $(.8\,\mathrm{PeV} \leq E < 2\,\mathrm{PeV})$")
    plot(sx, rad2deg(sy2), '^-',
         label="Theta (simulation, $E = 1\,\mathrm{PeV}$)")
    plot(x, rad2deg(y), '^-',
         label="Phi $(.8\,\mathrm{PeV} \leq E < 2\,\mathrm{PeV})$")
    plot(sx, rad2deg(sy), '^-',
         label="Phi (simulation, $E = 1\,\mathrm{PeV}$)")
    # Labels etc.
    xlabel("Core distance (m)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title(r"$\theta = 22.5 \pm %d^\circ, N_{MIP} = %d$" % (D2_Z, N))
    legend(loc='center right', numpoints=1)
    ylim(ymin=0)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-core-dist-phi-theta.pdf')

def plot_interarrival_times(h, k):
    f = lambda x, N, a: N * exp(a * x)

    figure()
    rcParams['text.usetex'] = False
    for shift in [-12, -13, -13.180220188, -14]:
        c = array(kascade_coincidences.search_coincidences(h, k, shift))
        n, bins, patches = hist(abs(c[:,0]) / 1e9, bins=linspace(0, 1, 200),
                                histtype='step', log=True,
                                label=r'$\Delta t = %.4f\,\mathrm{ns}$'
                                      % shift)

    b = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    popt, pcov = curve_fit(f, b, n)
    print "Interarrival times rate: %f" % popt[1]
    plot(b, f(b, *popt), label=r"$\lambda = %f$" % popt[1])

    xlabel("Time difference (s)")
    ylabel("Count")
    legend()
    ylim(ymin=10)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-interarrival-times.pdf')

    figure()
    rcParams['text.usetex'] = False
    shift = -13.180220188
    c = array(kascade_coincidences.search_coincidences(h, k[:20000], shift))
    l = len(c)
    n, bins, patches = hist(c[:,0] / 1e3, bins=linspace(-10, -5, 500),
                            histtype='step')
    xlabel("Time difference (us)")
    ylabel("Count")
    title(r"$\Delta t = %.9f\,\mathrm{s}$" % shift)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-interarrival-times-corr.pdf')



def negpos68(x):
    xmin = sorted(-x.compress(x <= 0))
    if len(xmin) == 0:
        return nan
    xmin = xmin[len(xmin) * 2 / 3]

    xmax = sorted(x.compress(x >= 0))
    if len(xmax) == 0:
        return nan
    xmax = xmax[len(xmax) * 2 / 3]

    return (xmin + xmax) / 2.

def get_raw_timings(events, i, j):
    return compress((events['n%d' % i] >= 1) & (events['n%d' % j] >= 1),
                    events['t%d' % i] - events['t%d' % j])

def plot_arrival_times_core(data, datasim):
    D2_Z = 5
    r = linspace(0, 100, 10)

    events = datasim.root.analysis.angle_0
    mr, tm, dt = [], [], []
    for r0, r1 in zip(r[:-1], r[1:]):
        sel = events.readWhere('(r0 <= r) & (r < r1)')
        t1 = sel['t1'] - sel['t2']
        t3 = sel['t3'] - sel['t2']
        t4 = sel['t4'] - sel['t2']
        t = array([t1.flatten(), t3.flatten(), t4.flatten()]).flatten()
        t = compress(-isnan(t), t)
        tm.append(mean(t))
        dt.append(negpos68(t))
        mr.append((r0 + r1) / 2.)
    mrsim, tmsim, dtsim = mr, tm, dt

    figure()
    rcParams['text.usetex'] = False
    events = data.root.reconstructions.raw
    mr, tm, dt = [], [], []
    for r0, r1 in zip(r[:-1], r[1:]):
        sel = events.readWhere('(r0 <= r) & (r < r1) & (%f <= k_theta) & '
                               '(k_theta < %f) & (%f <= k_energy) & '
                               '(k_energy < %f)' %
                               (-deg2rad(D2_Z), deg2rad(D2_Z), 8e14, 2e15))
        t = []
        for i, j in combinations((1, 2, 3, 4), 2):
            t += get_raw_timings(sel, i, j).tolist()
        t = array(t)

        tm.append(mean(t))
        dt.append(negpos68(t))
        mr.append((r0 + r1) / 2.)
    errorbar(mrsim, tmsim, yerr=dtsim, drawstyle='steps-mid', capsize=0,
             label="simulation")
    errorbar(mr, tm, yerr=dt, drawstyle='steps-mid', capsize=0,
             label="KASCADE")
    xlabel("Core distance (m)")
    ylabel("Mean arrival time (ns)")
    title(r"$\theta < %d^\circ, %.1f\,\mathrm{PeV} \leq E < "
          r"%.1f\,\mathrm{PeV}$" % (D2_Z, 8e14 / 1e15, 2e15 / 1e15))
    legend(numpoints=1, loc='best')
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-arrival-core-mean.pdf')

    figure()
    rcParams['text.usetex'] = False
    plot(mrsim, dtsim, 'v', label="simulation")
    plot(mr, dt, '^', label="KASCADE")
    xlabel("Core distance (m)")
    ylabel("arrival time spread (ns)")
    title(r"$\theta < %d^\circ, %.1f\,\mathrm{PeV} \leq E < "
          r"%.1f\,\mathrm{PeV}$" % (D2_Z, 8e14 / 1e15, 2e15 / 1e15))
    legend(numpoints=1, loc='best')
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-arrival-core-spread.pdf')

def plot_2d_results_phi(data):
    #table = data.getNode('/reconstructions', 'full_.9scaled_linear')
    table = data.root.reconstructions.full_linear
    events = table.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')

    mylog = vectorize(lambda x: log10(x) if x > 0 else 0)

    H, xedges, yedges = histogram2d(rad2deg(events['k_phi']),
                                    rad2deg(events['h_phi']), bins=200)
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2

    figure()
    rcParams['text.usetex'] = False
    gca().set_axis_bgcolor('b')
    contourf(x, y, H.T)
    colorbar()
    xlabel(r"$\phi_K\,(^\circ)$")
    ylabel(r"$\phi_H\,(^\circ)$")
    title("$N_{MIP} \geq 1$")
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-2d-phi.pdf')

    figure()
    rcParams['text.usetex'] = False
    gca().set_axis_bgcolor('b')
    lH = mylog(H)
    contourf(x, y, lH.T)
    colorbar()
    xlabel(r"$\phi_K\,(^\circ)$")
    ylabel(r"$\phi_H\,(^\circ)$")
    title("$N_{MIP} \geq 1$")
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-2d-phi-log.pdf')

    # Use the original timings (FSOT)
    #table = data.getNode('/reconstructions', 'full_.9scaled')
    table = data.root.reconstructions.full
    table2 = data.root.reconstructions.full_shifted
    events = table.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')
    events2 = table2.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')

    figure()
    rcParams['text.usetex'] = False
    phis = linspace(-pi, pi, 21)
    bins, mean_dphi, std_dphi, mean_dphi2, std_dphi2 = [], [], [], [], []
    for low, high in zip(phis[:-1], phis[1:]):
        # First table
        sel = events.compress((low <= events['k_phi']) &
                              (events['k_phi'] < high))
        dphi = sel['h_phi'] - sel['k_phi']
        dphi = (dphi + pi) % (2 * pi) - pi
        bins.append((low + high) / 2)
        mean_dphi.append(mean(dphi))
        std_dphi.append(std(dphi))

        # Second table
        sel = events2.compress((low <= events2['k_phi']) &
                              (events2['k_phi'] < high))
        dphi = sel['h_phi'] - sel['k_phi']
        dphi = (dphi + pi) % (2 * pi) - pi
        mean_dphi2.append(mean(dphi))
        std_dphi2.append(std(dphi))
    errorbar(rad2deg(bins), rad2deg(mean_dphi), yerr=rad2deg(std_dphi),
             #drawstyle='steps-mid', fmt='-o', label='Uncorrected')
             fmt='o', label='Uncorrected')
    errorbar(rad2deg(bins), rad2deg(mean_dphi2), yerr=rad2deg(std_dphi2),
             #drawstyle='steps-mid', fmt='-o', label='Time-shifted')
             fmt='o', label='Time-shifted')
    axhline(0, c='k')
    xlabel(r"$\phi_K\,(^\circ)$")
    ylabel(r"Mean $\phi_H - \phi_K\,(^\circ)$")
    title("$N_{MIP} \geq 1$")
    legend(numpoints=1, loc='best')
    xlim(-180, 180)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-bin-phi.pdf')

def plot_2d_results_theta(data):
    #table = data.getNode('/reconstructions', 'full_.9scaled_linear')
    table = data.root.reconstructions.full_linear
    events = table.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')

    H, xedges, yedges = histogram2d(rad2deg(events['k_theta']),
                                    rad2deg(events['h_theta']),
                                    bins=[linspace(0, 40, 201),
                                          linspace(0, 40, 201)])

    figure()
    gca().set_axis_bgcolor('b')
    rcParams['text.usetex'] = False
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2
    contourf(x, y, H.T)
    colorbar()
    xlabel(r"$\theta_K\,(^\circ)$")
    ylabel(r"$\theta_H\,(^\circ)$")
    title("$N_{MIP} \geq 1$")
    axis('auto')
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-2d-theta.pdf')

    # Use the original (FSOT) timings
    #table = data.getNode('/reconstructions', 'full_.9scaled')
    table = data.root.reconstructions.full
    table2 = data.root.reconstructions.full_shifted
    events = table.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')
    events2 = table2.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')

    figure()
    rcParams['text.usetex'] = False
    thetas = linspace(0, deg2rad(40), 21)
    bins, mean_dtheta, std_dtheta, mean_dtheta2, std_dtheta2 = [], [], \
                                                               [], [], []
    for low, high in zip(thetas[:-1], thetas[1:]):
        # First table
        sel = events.compress((low <= events['k_theta']) &
                              (events['k_theta'] < high))
        dtheta = sel['h_theta'] - sel['k_theta']
        bins.append((low + high) / 2)
        mean_dtheta.append(mean(dtheta))
        std_dtheta.append(std(dtheta))

        # Second table
        sel = events2.compress((low <= events2['k_theta']) &
                              (events2['k_theta'] < high))
        dtheta = sel['h_theta'] - sel['k_theta']
        mean_dtheta2.append(mean(dtheta))
        std_dtheta2.append(std(dtheta))
    errorbar(rad2deg(bins), rad2deg(mean_dtheta), yerr=rad2deg(std_dtheta),
             #drawstyle='steps-mid', fmt='-o', label="Uncorrected")
             fmt='o', label="Uncorrected")
    errorbar(rad2deg(bins), rad2deg(mean_dtheta2),
             yerr=rad2deg(std_dtheta2),
             #drawstyle='steps-mid', fmt='-o', label="Time-shifted")
             fmt='o', label="Time-shifted")
    axhline(0, c='k')
    xlabel(r"$\theta_K\,(^\circ)$")
    ylabel(r"Mean $\theta_H - \theta_K\,(^\circ)$")
    title("$N_{MIP} \geq 1$")
    legend(loc='best', numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-bin-theta.pdf')

def plot_energy_zenith_bin(data, tablename):
    table = data.getNode('/reconstructions', tablename)


    figure()
    rcParams['text.usetex'] = False
   
    energies = arange(14, 17.1, 1.)
    thetas = linspace(0, deg2rad(40), 21)

    for elow, ehigh in zip(energies[:-1], energies[1:]):
        events = table.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1) & '
                                 '(%f <= k_energy) & (k_energy < %f)' %
                                 (10 ** elow, 10 ** ehigh))
        print elow, ehigh, len(events)
        bins, mean_dtheta, std_dtheta, mean_dtheta2, std_dtheta2 = \
            [], [], [], [], []
        for low, high in zip(thetas[:-1], thetas[1:]):
            # First table
            sel = events.compress((low <= events['k_theta']) &
                                  (events['k_theta'] < high))
            dtheta = sel['h_theta'] - sel['k_theta']
            bins.append((low + high) / 2)
            mean_dtheta.append(mean(dtheta))
            std_dtheta.append(std(dtheta))

        plot(rad2deg(bins), rad2deg(mean_dtheta), 'o',
             label="%.1f <= log(E) < %.1f" % (elow, ehigh))
    axhline(0, c='k')
    xlabel(r"$\theta_K\,(^\circ)$")
    ylabel(r"Mean $\theta_H - \theta_K\,(^\circ)$")
    title("$N_{MIP} \geq 1$")
    legend(loc='best', numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-bin-theta-energy.pdf')

def plot_zenith_bin_12(data):
    table = data.root.reconstructions.full3_linear
    events = table.readWhere('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')

    figure()
    rcParams['text.usetex'] = False
    thetas = linspace(0, deg2rad(40), 21)
    bins, mean_dtheta, mean_dtheta2 = [], [], []
    for low, high in zip(thetas[:-1], thetas[1:]):
        sel = events.compress((low <= events['k_theta']) &
                              (events['k_theta'] < high))
        dtheta = sel['h_theta1'] - sel['k_theta']
        dtheta2 = sel['h_theta2'] - sel['k_theta']
        print (low + high) / 2, mean(dtheta), mean(dtheta2)
        bins.append((low + high) / 2)
        mean_dtheta.append(mean(dtheta))
        mean_dtheta2.append(mean(dtheta2))

    plot(rad2deg(bins), rad2deg(mean_dtheta), 'o', label="Theta 1")
    plot(rad2deg(bins), rad2deg(mean_dtheta2), 'o', label="Theta 2")
    axhline(0, c='k')
    xlabel(r"$\theta_K\,(^\circ)$")
    ylabel(r"Mean $\theta_H - \theta_K\,(^\circ)$")
    title("$N_{MIP} \geq 1$")
    legend(loc='best', numpoints=1)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-bin-theta-12.pdf')

def reconstruct_optimal_coincidences(h, k, initial, start=None, limit=None):
    """Determine optimal timeshift for coincidences and reconstruct"""

    t0 = max(h[0][1], k[0][1])

    if start: 
        start *= int(1e9)
        for i, t in enumerate([x[1] for x in h]):
            if t - t0 > start:
                break
        h = h[i:]
        for i, t in enumerate([x[1] for x in k]):
            if t - t0 > start:
                break
        k = k[i:]
    
    t0 += start

    if limit:
        limit *= int(1e9)
        for i, t in enumerate([x[1] for x in h]):
            if t - t0 > limit:
                break
        h = h[:i]
        for i, t in enumerate([x[1] for x in k]):
            if t - t0 > limit:
                break
        k = k[:i]

    #h = [(x[0], x[1] ^ (2 ** 10 - 1)) for x in h]
    #h = [(x[0], x[1] ^ (2 ** 11 - 1)) for x in h]
    #c = kascade_coincidences.search_coincidences(h, k, initial)
    #m = median([x[0] for x in c])
    #shift = initial - m / 1e9
    #print shift

    shift = initial
    c = kascade_coincidences.search_coincidences(h, k, shift)
    m = median([x[0] for x in c])
    print "Median of residual time differences (ns):", m
    print "Length: ", len(c)
    #figure()
    hist([x[0] for x in c], bins=range(-2000, 5000, 5), histtype='step')

    # Drop outliers
    c = [x for x in c if -1e6 < x[0] < 1e6]

    #figure()
    rcParams['text.usetex'] = False
    #plot([(k[x[2]][1] - t0) / 1e9 for x in c], [x[0] for x in c], ',')

    dt = [x[0] for x in c]
    h_idx = [x[1] for x in c]
    k_idx = [x[2] for x in c]
    h_t = [h[x][1] % int(1e9) for x in h_idx]
    k_t = [k[x][1] % int(1e9) for x in k_idx]

    #plot(h_t, dt, 'o-')

    xlabel("Time (s)")
    ylabel("Residual time difference (ns)")
    axis('tight')

def time_plot(h, k, initial, batchsize=5000, limit=None):
    tl = []
    ml = []

    t0 = max(h[0][1], k[0][1])
    if limit:
        limit *= int(1e9)
        for i, t in enumerate([x[1] for x in h]):
            if t - t0 > limit:
                break
        h = h[:i]
        for i, t in enumerate([x[1] for x in k]):
            if t - t0 > limit:
                break
        k = k[:i]

    c = kascade_coincidences.search_coincidences(h, k, initial, dtlimit=100000)

    figure()
    for i in range(0, len(c), batchsize):
        tc = c[i:i + batchsize]
        
        print array([x[0] for x in tc])
        m = median([x[0] for x in tc])
        #t = [x[0] - m for x in tc]
        t = [x[0] for x in tc]
        ml.append(m)
        tl.append(h[tc[0][1]][1] - t0)
        print "Median:", m
        print "Duration: %.1f hours" % ((h[tc[-1][1]][1] - h[tc[0][1]][1]) / 1e9 / 3600)
        print "Length:", len(t)

        hist(t, bins=range(-2000, 5000, 5), histtype='step')

    xlabel("Time (ns)")
    ylabel("Count")
    axis('tight')

    figure()
    plot([x / 1e9 / 86400 for x in tl], ml, '-')
    xlabel("Time (days)")
    ylabel("Median (ns)")


if __name__ == '__main__':
    # invalid values in arcsin will be ignored (nan handles the situation
    # quite well)
    np.seterr(invalid='ignore')

    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'a')

    try:
        h, k
    except NameError:
        try:
            print "Reading timing datasets from disk"
            h = data.root.datasets.h.read()
            k = data.root.datasets.k.read()
        except tables.NoSuchNodeError:
            print "Building timing datasets from disk"
            h, k = kascade_coincidences.get_arrays_from_tables(
                        data.root.hisparc.cluster_kascade.station_601.events,
                        data.root.kascade.events)
            data.createGroup('/', 'datasets')
            data.createArray('/datasets', 'h', h)
            data.createArray('/datasets', 'k', k)

    if 'coincidences' in data.root:
        print "Reading coincidences from disk"
        events = data.root.coincidences.events
    else:
        print "Searching for coincidences"
        c = kascade_coincidences.search_coincidences(h, k, -13.180220188,
                                                     dtlimit=1e9)
        print "Storing data on disk"
        data.createGroup('/', 'coincidences',
                         "HiSPARC / KASCADE coincidences")
        events = data.createTable('/coincidences', 'events', Coincidence)
        kascade_coincidences.store_coincidences(events,
                                data.root.hisparc.cluster_kascade.station_601.events,
                                data.root.kascade.events, c)

    try:
        timing_data
        timing_data_linear
    except NameError:
        try:
            print "Reading timing data from disk"
            timing_data = data.root.analysis.timing_data.read()
            timing_data_linear = data.root.analysis.timing_data_linear.read()
        except tables.NoSuchNodeError:
            print "Analysing timing data"
            timing_data, timing_data_linear = process_traces(events,
                        data.root.hisparc.cluster_kascade.station_601.blobs)
            print "Writing timing data to disk"
            data.createGroup('/', 'analysis')
            data.createArray('/analysis', 'timing_data', timing_data)
            data.createArray('/analysis', 'timing_data_linear',
                             timing_data_linear)

    try:
        sim_data
    except NameError:
        sim_data = tables.openFile('data-e15.h5', 'r')


    #print "Reconstructing events..."
    #reconstruct_angles(data, 'full', events, timing_data)
    #reconstruct_angles(data, 'full_linear', events, timing_data_linear)
    #reconstruct_angles(data, 'full_shifted', events, timing_data,
    #                   [0.25, 0, 1.17, -0.21])
    #reconstruct_angles(data, 'full_shifted_linear', events,
    #                   timing_data_linear, [0.26, 0, 1.20, -0.19])

    #do_reconstruction_plots(data, 'full', 'full_linear', sim_data, 'full')
    #plot_interarrival_times(h, k)
    
    #time_plot(h, k, initial=-13.180212844, batchsize=5000, limit=10 * 86400)
    #time_plot(h, k, initial=-13.180212844, batchsize=5000, limit=86400)

    time_plot(h, k, initial=-15, batchsize=5000, limit=None)
