import tables
import csv
import zlib

from pylab import *
from scipy.optimize import curve_fit

from scipy import integrate
from scipy.special import erf

import progressbar

DATAFILE = 'kascade.h5'

USE_TEX = False


D_Z = 1 # Delta zenith (used in data selection)

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


def do_reconstruction_plots(data, tablename, table2name, sim_data,
                            sim_tablename):
    """Make plots based upon earlier reconstructions"""

    table = data.getNode('/reconstructions', tablename)
    table2 = data.getNode('/reconstructions', table2name)
    sim_table = sim_data.getNode('/reconstructions', sim_tablename)

    #plot_2d_results_phi(data)
    #plot_2d_results_theta(data)

    #plot_uncertainty_mip(table, sim_table)
    plot_uncertainty_zenith(table, sim_table)
    #plot_uncertainty_zenith2(table, table2)
    #plot_uncertainty_energy(table)
    #plot_mip_core_dists_mean(table, sim_table)
    #plot_zenith_core_dists_mean(table, sim_table)
    #plot_uncertainty_core_dist_phi_theta(table, sim_table)

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
    #savefig('plots/auto-results-MIP.pdf')
    print

    with open('plotdata/auto-results-MIP_marks.dat', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(('mip', 'theta', 'theta_sim', 'phi', 'phi_sim'))
        writer.writerows(zip(x, rad2deg(y), rad2deg(sy), rad2deg(y2),
                             rad2deg(sy2)))
    with open('plotdata/auto-results-MIP_lines.dat', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(('mip', 'theta', 'phi'))
        writer.writerows(zip(cx, rad2deg(cy2), rad2deg(cy)))

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

    with open('plotdata/auto-results-zenith_marks_theta.dat', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(('zenith', 'theta', 'theta_sim'))
        writer.writerows(zip(rad2deg(x), rad2deg(y), rad2deg(sy)))
    with open('plotdata/auto-results-zenith_marks_phi.dat', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(('zenith', 'phi', 'phi_sim'))
        writer.writerows(zip(rad2deg(x[1:]), rad2deg(y2[1:]),
                             rad2deg(sy2[1:])))
    with open('plotdata/auto-results-zenith_lines.dat', 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(('zenith', 'theta', 'phi', 'phi_sintheta'))
        writer.writerows(zip(rad2deg(cx), rad2deg(cy2), rad2deg(cy),
                             rad2deg(cy3)))

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

    bins = linspace(0, 1, 200)
    ns = []
    for shift in [-12, -13, -13.180220188, -14]:
        c = array(kascade_coincidences.search_coincidences(h, k, shift))
        n, bins = histogram(abs(c[:, 0]) / 1e9, bins=bins)
        n = n.tolist() + [n[-1]]
        n = [u if u else .1 for u in n]
        ns.append(n)
    #with open('plotdata/interarrival_times.dat', 'w') as f:
    #    writer = csv.writer(f, delimiter='\t')
    #    writer.writerow(('bin', 'n12', 'n13', 'n13x', 'n14'))
    #    writer.writerows(zip(bins, *ns))

    figure()
    rcParams['text.usetex'] = False
    for shift in [-12, -13, -13.180220188, -14]:
        c = array(kascade_coincidences.search_coincidences(h, k, shift))
        n, bins, patches = hist(abs(c[:, 0]) / 1e9, bins=linspace(0, 1, 200),
                                histtype='step', log=True,
                                label=r'$\Delta t = %.4f\,\mathrm{ns}$'
                                      % shift)

    b = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    popt, pcov = curve_fit(f, b, n)
    print "Interarrival times rate: %f" % popt[1]
    print "Scaling factor: %f" % popt[0]
    plot(b, f(b, *popt), label=r"$\lambda = %f$" % popt[1])

    xlabel("Time difference (s)")
    ylabel("Count")
    legend()
    ylim(ymin=10)
    if USE_TEX:
        rcParams['text.usetex'] = True
    savefig('plots/auto-results-interarrival-times.pdf')

#    figure()
#    rcParams['text.usetex'] = False
#    shift = -13.180220188
#    c = array(kascade_coincidences.search_coincidences(h, k[:20000], shift))
#    l = len(c)
#    n, bins, patches = hist(c[:,0] / 1e3, bins=linspace(-10, -5, 500),
#                            histtype='step')
#    xlabel("Time difference (us)")
#    ylabel("Count")
#    title(r"$\Delta t = %.9f\,\mathrm{s}$" % shift)
#    if USE_TEX:
#        rcParams['text.usetex'] = True
#    savefig('plots/auto-results-interarrival-times-corr.pdf')



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

    global H, x, y
    H, xedges, yedges = histogram2d(rad2deg(events['k_phi']),
                                    rad2deg(events['h_phi']), bins=20)
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2

    # plotdata output
    m = H.max()
    wx = x[1] - x[0]
    wy = y[1] - y[0]
    with open('plotdata/auto-results-2d-phi.inp', 'w') as f:
        for cx, u in zip(x, H):
            for cy, v in zip(y, u):
                if v != 0:
                    nv = v / m
                    rx = wx * nv
                    ry = wy * nv
                    x0 = cx - .5 * rx
                    x1 = cx + .5 * rx
                    y0 = cy - .5 * ry
                    y1 = cy + .5 * ry
                    f.write('\\path[data] (axis cs: %f, %f) rectangle '
                            '(axis cs: %f, %f);\n' % (x0, y0, x1, y1))
    return

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
                                    bins=[linspace(0, 40, 21),
                                          linspace(0, 40, 21)])
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2

    # plotdata output
    m = H.max()
    wx = x[1] - x[0]
    wy = y[1] - y[0]
    with open('plotdata/auto-results-2d-theta.inp', 'w') as f:
        for cx, u in zip(x, H):
            for cy, v in zip(y, u):
                if v != 0:
                    nv = v / m
                    rx = wx * nv
                    ry = wy * nv
                    x0 = cx - .5 * rx
                    x1 = cx + .5 * rx
                    y0 = cy - .5 * ry
                    y1 = cy + .5 * ry
                    f.write('\\path[data] (axis cs: %f, %f) rectangle '
                            '(axis cs: %f, %f);\n' % (x0, y0, x1, y1))
    return

    figure()
    gca().set_axis_bgcolor('b')
    rcParams['text.usetex'] = False
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
