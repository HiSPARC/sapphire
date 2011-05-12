from __future__ import division

import tables
import progressbar as pb

from build_dataset import DETECTORS
from plot_utils import *

from pylab import *
from scipy.optimize import curve_fit, fmin
from scipy import stats

from analysis_kascade import std_t, calc_phi

from landau import Scintillator, discrete_convolution


DATAFILE = 'kascade.h5'
GROUP = '/efficiency'


phi1 = calc_phi(1, 3)
phi2 = calc_phi(1, 4)


def progressbar(*args, **kwargs):
    progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])
    return progress(*args, **kwargs)

def main():
    #plot_pulseheight_trigger_histogram()

    #plot_density_histogram()
    #plot_energy_histogram()
    #plot_core_distance()
    #plot_core_alpha()
    #plot_core_pos()
    #plot_zenith()
    #plot_azimuth()
    #plot_T200()
    #plot_P200()
    #plot_particle_density()
    #plot_density_fraction_err()
    #plot_density_fraction()
    #optimize_trigger_prob()
    #plot_theta_phi_corr()

    #plot_density_spread()
    #plot_density_coredist()

    #kevents, hevents, labels, sfx = get_regions_data()
    #plot_regions_core_pos(kevents, hevents, labels, sfx)
    #plot_regions_energy_count(kevents, hevents, labels, sfx)
    #plot_regions_energy_ratio(kevents, hevents, labels, sfx)
    #plot_regions_zenith(kevents, hevents, labels, sfx)
    #plot_regions_azimuth(kevents, hevents, labels, sfx)
    #plot_regions_T200(kevents, hevents, labels, sfx)
    #plot_regions_P200(kevents, hevents, labels, sfx)

    #kevents, hevents, labels, sfx = get_p_regions_data()
    #plot_regions_energy_count(kevents, hevents, labels, sfx)
    #plot_regions_energy_ratio(kevents, hevents, labels, sfx)
    #plot_regions_P200(kevents, hevents, labels, sfx)
    #plot_regions_core_pos(kevents, hevents, labels, sfx)

    #plot_reconstruction_efficiency_zenith()
    #plot_reconstruction_efficiency_azimuth()
    #plot_reconstruction_efficiency_density()
    #plot_reconstruction_efficiency_coredist()

    #plot_detector_efficiency()

    #plot_rec_eff_timings()
    #plot_rvalue_core()

    #electron_gamma_trigger()

    #plot_poisson_trigger_cuts()

    #plot_charged_particles_poisson()
    plot_conv_poisson()

    pass

def plot_pulseheight_trigger_histogram():
    events = data.root.hisparc.cluster_kascade.station_601.events

    figure()
    hist(events.col('pulseheights') * .57, bins=linspace(0, 1600, 201),
         log=True, histtype='step')
    xlabel("Pulseheight (mV)")
    ylabel("Count")
    ylim(ymin=10)
    legend(['Detector %d' % x for x in range(1, 5)])
    axvline(30)
    axvline(70)
    savefig("plots/pulseheight_trigger_hist.pdf")

def plot_density_histogram():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['k_dens_e'], bins=linspace(0, 10, 200),
         histtype='step')
    xlabel("Electron density (m$^{-1}$)")
    ylabel("Count")
    savefig("plots/density_histogram.pdf")

def plot_energy_histogram():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['k_energy'] / 1e15, bins=logspace(-1, 2, 200),
         histtype='step', label="KASCADE trigger", normed=True)
    hist(events.readWhere('self_triggered == True')['k_energy'] / 1e15,
         bins=logspace(-1, 2, 200), histtype='step',
         label="HiSPARC trigger", normed=True)
    xlabel("Energy (PeV)")
    xscale('log')
    ylabel("Count")
    legend(loc='best')
    savefig("plots/energy_histogram.pdf")

def plot_core_distance():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['core_dist'], bins=linspace(0, 200, 200),
         histtype='step', label="KASCADE trigger", normed=True)
    hist(events.readWhere('self_triggered == True')['core_dist'],
         bins=linspace(0, 200, 200), histtype='step',
         label="HiSPARC trigger", normed=True)
    axvline(15)
    axvline(30)
    xlabel("Core distance (m)")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/core_distance.pdf")
    
    figure()
    hk, bins = histogram(events[:]['core_dist'], bins=linspace(0, 200,
                                                               200))
    hh, bins = histogram(
            events.readWhere('self_triggered == True')['core_dist'],
            bins=linspace(0, 200, 200))
    x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
    hr = 1. * hh / hk
    plot(x, hr)
    axvline(15)
    axvline(30)
    xlabel("Core distance (m)")
    ylabel("Ratio of selftriggers")
    savefig("plots/core_distance_ratio.pdf")

def plot_core_alpha():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['core_alpha'] / pi, bins=200, histtype='step',
         label="KASCADE trigger", normed=True)
    hist(events.readWhere('self_triggered == True')['core_alpha'] / pi,
         bins=200, histtype='step', label="HiSPARC trigger", normed=True)
    xlabel("Ground angle to core ($\mathrm{rad} \pi^{-1}$)")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/core_alpha.pdf")

def plot_core_pos():
    events = data.getNode(GROUP, 'events')

    figure()
    cores = events[:]['k_core_pos']
    contour_histogram2d(cores[:,0], cores[:,1], bins=200)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions")
    savefig("plots/core_pos.pdf")

    figure()
    contour_histogram2d(cores[:,0], cores[:,1], bins=200, log=True)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions (log)")
    savefig("plots/core_pos_log.pdf")
    
    figure()
    cores = events.readWhere('self_triggered == True')['k_core_pos']
    contour_histogram2d(cores[:,0], cores[:,1], bins=200)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions (self-triggered)")
    savefig("plots/core_pos_trig.pdf")

    figure()
    contour_histogram2d(cores[:,0], cores[:,1], bins=200, log=True)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions (self-triggered) (log)")
    savefig("plots/core_pos_trig_log.pdf")

def plot_zenith():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(rad2deg(events[:]['k_zenith']), bins=200, histtype='step',
         label="KASCADE trigger", normed=True)
    hist(rad2deg(events.readWhere('self_triggered == True')['k_zenith']),
         bins=200, histtype='step', label="HiSPARC trigger", normed=True)
    xlabel("Zenith $(^\circ)$")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/zenith_histogram.pdf")
   
    figure()
    hk, bins = histogram(rad2deg(events[:]['k_zenith']), bins=200)
    hh, bins = histogram(
            rad2deg(events.readWhere('self_triggered == True')['k_zenith']),
            bins=200)
    x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
    hr = 1. * hh / hk
    plot(x, hr)
    xlabel("Zenith $(^\circ)$")
    ylabel("Ratio of self-triggers")
    legend(loc='best')
    savefig("plots/zenith_ratio.pdf")

def plot_azimuth():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(rad2deg(events[:]['k_azimuth']), bins=200, histtype='step',
         label="KASCADE trigger", normed=True)
    hist(rad2deg(events.readWhere('self_triggered == True')['k_azimuth']),
         bins=200, histtype='step', label="HiSPARC trigger", normed=True)
    xlabel("Azimuth $(^\circ)$")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/azimuth_histogram.pdf")
   
    figure()
    hk, bins = histogram(rad2deg(events[:]['k_azimuth']), bins=200)
    hh, bins = histogram(
            rad2deg(events.readWhere('self_triggered == True')['k_azimuth']),
            bins=200)
    x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
    hr = 1. * hh / hk
    plot(x, hr)
    xlabel("Azimuth $(^\circ)$")
    ylabel("Ratio of self-triggers")
    legend(loc='best')
    savefig("plots/azimuth_ratio.pdf")

def plot_T200():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['k_T200'], bins=100, histtype='step',
         label="KASCADE trigger", normed=True)
    hist(events.readWhere('self_triggered == True')['k_T200'], bins=200,
         histtype='step', label="HiSPARC trigger", normed=True)
    xlabel("Temperature $(^\circ\!\mathrm{C})$")
    ylabel("Count")
    # FIXME
    ylim(ymax=.15)
    legend(loc='best')
    savefig("plots/T200_histogram.pdf")

    figure()
    hk, bins = histogram(events[:]['k_T200'], bins=200)
    hh, bins = histogram(
            events.readWhere('self_triggered == True')['k_T200'],
            bins=200)
    x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
    hr = 1. * hh / hk
    plot(x, hr)
    xlabel("Temperature $(^\circ\!\mathrm{C})$")
    ylabel("Ratio of self-triggers")
    legend(loc='best')
    savefig("plots/T200_ratio.pdf")

def plot_P200():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['k_P200'], bins=200, histtype='step',
         label="KASCADE trigger", normed=True)
    hist(events.readWhere('self_triggered == True')['k_P200'], bins=200,
         histtype='step', label="HiSPARC trigger", normed=True)
    xlabel("Pressure (hPa)")
    ylabel("Count")
    # FIXME
    ylim(ymax=.2)
    legend(loc='best')
    savefig("plots/P200_histogram.pdf")
   
    figure()
    hk, bins = histogram(events[:]['k_P200'], bins=200)
    hh, bins = histogram(
            events.readWhere('self_triggered == True')['k_P200'],
            bins=200)
    x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
    hr = 1. * hh / hk
    plot(x, hr)
    xlabel("Pressure (hPa)")
    ylabel("Ratio of self-triggers")
    legend(loc='best')
    savefig("plots/P200_ratio.pdf")

def plot_particle_density():
    events = data.getNode(GROUP, 'events')

    figure()
    edges = linspace(0, 200, 200)
    x, y, y2 = [], [], []
    for r1, r2 in zip(edges[:-1], edges[1:]):
        e = events.readWhere('(r1 <= core_dist) & (core_dist < r2)')
        x.append((r1 + r2) / 2)
        y.append(median(e[:]['k_dens_e']))
        y2.append(median(e[:]['k_dens_mu']))
    plot(x, y, label="Electrons")
    plot(x, y2, label="Muons")
    xlabel("Core distance (m)")
    ylabel("Median electron density $(\mathrm{m}^{-2})$")
    legend(loc='best')
    savefig("plots/core_dist_pdens.pdf")

def plot_density_fraction_err():
    global kevents
    events = data.getNode(GROUP, 'events')

    #q = '(k_Num_e > 10 ** 5)'
    q = '(k_Num_e > 10 ** 4.6) & (k_Num_e < 10 ** 7) & (core_dist_center <= 90) & (k_theta <= %f)' % deg2rad(30)

    clf()
    bins = linspace(0, 10, 201)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    kevents = events.readWhere(q)
    hevents = events.readWhere(q + '& (self_triggered == True)')
    kdens = kevents['k_dens_e']
    kdens_m = mean(kdens, 1)
    hdens = hevents['k_dens_e']
    hdens_m = mean(hdens, 1)

    err = [[] for u in xrange(len(bins) - 1)]
    for idx, v in zip(bins.searchsorted(hdens_m), hdens):
        try:
            err[idx - 1].extend(v.tolist())
        except IndexError:
            pass
    errb = [std(u) for u in err]

    hk, bins = histogram(kdens_m, bins=bins)
    hh, bins = histogram(hdens_m, bins=bins)
    hr = 1. * hh / hk
    errorbar(x, hr, xerr=errb, label="Data", capsize=0)
    
    # Poisson probability of zero particles in detector
    p0 = exp(-.5 * x)
    p = 1 - (4 * (1 - p0) * p0 ** 3 + p0 ** 4)
    plot(x, p, label="Poisson")

    xlabel("Electron density (m$^{-1}$)")
    ylabel("Probability of trigger")
    legend(loc='best')
    xlim(0, 10)
    ylim(0, 1.05)
    savefig("plots/density_fraction_err.pdf")

def plot_density_fraction():
    global events
    events = data.getNode(GROUP, 'events')

    q = '(k_Num_e > 1e5)'

    figure()
    bins = linspace(0, 10, 201)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    for Z in ['k_zenith < %.2f' % deg2rad(40), 'k_zenith < %.2f' %
              deg2rad(5), 'k_zenith > %.2f' % deg2rad(25)]:
        kevents = events.readWhere(q + '& (%s)' % Z)
        print Z, mean(log10(kevents['k_energy']))
        hevents = events.readWhere(q + '& (%s) & (self_triggered == True)' % Z)
        #hevents = events.readWhere('(%s) & (n_high >= 3)' % Z)
        kdens = mean(kevents['k_dens_e'] + kevents['k_dens_mu'], 1)
        hdens = mean(hevents['k_dens_e'] + hevents['k_dens_mu'], 1)

        subplot(211)
        hk, bins = histogram(kdens, bins=bins)
        hh, bins = histogram(hdens, bins=bins)
        hr = 1. * hh / hk
        plot(x, hr, label=Z)

        subplot(212)
        kdens *= cos(kevents['k_zenith'])
        hdens *= cos(hevents['k_zenith'])
        hk, bins = histogram(kdens, bins=bins)
        hh, bins = histogram(hdens, bins=bins)
        hr = 1. * hh / hk
        plot(x, hr, label=Z)

    # Poisson probability of zero particles in detector
    p0 = exp(-.5 * x)
    p = 1 - (4 * (1 - p0) * p0 ** 3 + p0 ** 4)
    #p = 1 - (4 * (1 - p0) * p0 ** 3 + p0 ** 4 + 6 * (1 - p0) ** 2 * p0 ** 2)

    subplot(211)
    plot(x, p, label="Poisson")
    ylabel("Probability of trigger")
    legend(loc='best')
    title("Without correction")

    subplot(212)
    plot(x, p, label="Poisson")
    ylabel("Probability of trigger")
    legend(loc='best')
    xlabel("Electron density (m$^{-1}$)")
    title(r"$\cos\theta$ correction")
    savefig("plots/density_trigger_probability.pdf")

def optimize_trigger_prob():
    events = data.getNode(GROUP, 'events')

    figure()
    bins = linspace(0, 10, 201)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    hrs = []
    for Z in ['k_zenith < %.2f' % deg2rad(40), 'k_zenith < %.2f' %
              deg2rad(10), 'k_zenith > %.2f' % deg2rad(25)]:
        kevents = events.readWhere(Z)
        hevents = events.readWhere('(%s) & (self_triggered == True)' % Z)
        #hevents = events.readWhere('(%s) & (n_high >= 3)' % Z)
        kdens = median(kevents['k_dens_e'], 1) * cos(kevents['k_zenith'])
        hdens = median(hevents['k_dens_e'], 1) * cos(hevents['k_zenith'])

        hk, bins = histogram(kdens, bins=bins)
        hh, bins = histogram(hdens, bins=bins)
        hr = 1. * hh / hk
        hrs.append(hr)

    # Poisson probability of zero particles in detector
    p = lambda x, c: 1 - (4 * (1 - exp(-.5 * (x + c))) * \
                     exp(-.5 * (x + c)) ** 3 + exp(-.5 * (x + c)) ** 4)
                     #+ 6 * (1 - exp(-.5 * (x + c))) ** 2 * \
                     #exp(-.5 * (x + c)) ** 2)

    popt, pcov = curve_fit(p, x, hrs[0])
    for y in hrs:
        chisq = sum([1. / (.5 * u) * (v - w) ** 2 for u, v, w in
                    zip(x, y, p(x, popt[0])) if u > .25])
        plot(x + popt[0], y, label='chisq: %.2e' % chisq)

    plot(x, p(x, 0), label="Poisson")
    xlabel("Electron density (m$^{-1}$)")
    ylabel("Probability of trigger")
    ylim(ymax=1.05)
    title(r"$\cos\theta$ correction, shift: %.2f" % popt[0])
    legend(loc='best')
    savefig("plots/optimize_trigger_prob.pdf")

def plot_theta_phi_corr():
    events = data.root.reconstructions.full_linear[:]
    events = events.compress(events['k_dens_e'][:,1] < .3)

    figure()
    plot(rad2deg(events['k_phi']), rad2deg(events['h_phi']), ',')
    xlabel("KASCADE phi ($^\circ$)")
    ylabel("HiSPARC phi ($^\circ$)")
    savefig("plots/low_dens_phi_corr.pdf")

    figure()
    plot(rad2deg(events['k_theta']), rad2deg(events['h_theta']), ',')
    xlabel("KASCADE theta ($^\circ$)")
    ylabel("HiSPARC theta ($^\circ$)")
    savefig("plots/low_dens_theta_corr.pdf")

def plot_density_spread():
    events = data.getNode(GROUP, 'events')
    events = events.readWhere('(.9e15 <= k_energy) & (k_energy < 1.1e15)')

    dmean = array([mean(u['k_dens_e']) for u in events[:]])
    dspread = array([max(u['k_dens_e']) - min(u['k_dens_e']) for u in
                     events[:]])

    bins = linspace(0, 10, 101)
    xspread = [compress((u <= dmean) & (dmean < v), dspread) for
                     u, v in zip(bins[:-1], bins[1:])]
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    y = array([mean(u) for u in xspread])
    yerr = array([std(u) for u in xspread])

    figure()
    errorbar(x, y, yerr=yerr, capsize=0)
    xlabel("Mean density (m$^{-1}$)")
    ylabel("Density spread (m$^{-1}$)")

def plot_density_coredist():
    events = data.getNode(GROUP, 'events')
    events = events.readWhere('(.9e15 <= k_energy) & (k_energy < 1.1e15)'
                              ' & (k_zenith < .35)')

    dmean = array([mean(u['k_dens_e']) for u in events[:]])
    coredist = events[:]['core_dist']

    bins = linspace(0, 10, 101)
    xcoredist = [compress((u <= dmean) & (dmean < v), coredist) for u, v
                 in zip(bins[:-1], bins[1:])]
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    y = array([mean(u) for u in xcoredist])
    yerr = array([std(u) for u in xcoredist])

    figure()
    errorbar(x, y, yerr=yerr, capsize=0)
    xlabel("Mean density (m$^{-1}$)")
    ylabel("Mean core distance (m)")

def plot_regions_core_pos(kevents, hevents, labels, sfx):
    for i, (k, h, label) in enumerate(zip(kevents, hevents, labels)):
        figure()
        cores = k[:]['k_core_pos']
        contour_histogram2d(cores[:,0], cores[:,1], bins=200)
        plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
        xlabel("(m)")
        ylabel("(m)")
        title("Shower core positions (%s)" % label)
        savefig("plots/region_core_pos_%d_%s.pdf" % (i, sfx))

    for i, (k, h, label) in enumerate(zip(kevents, hevents, labels)):
        figure()
        cores = h[:]['k_core_pos']
        contour_histogram2d(cores[:,0], cores[:,1], bins=200)
        plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
        xlabel("(m)")
        ylabel("(m)")
        title("Shower core positions (self-triggered) (%s)" % label)
        savefig("plots/region_core_pos_trig_%d_%s.pdf" % (i, sfx))

def plot_regions_energy_count(kevents, hevents, labels, sfx):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(k['k_energy'] / 1e15, bins=logspace(-1, 2, 200),
             histtype='step', label=label, normed=True)
        subplot(212)
        hist(h['k_energy'] / 1e15, bins=logspace(-1, 2, 200),
             histtype='step', label=label, normed=True)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    xscale('log')
    ylabel("Count")
    legend(loc='best')
    subplot(212)
    xlabel("Energy (PeV)")
    xscale('log')
    ylabel("Count")
    legend(loc='best')
    savefig("plots/region_energy_%s.pdf" % sfx)

def plot_regions_energy_ratio(kevents, hevents, labels, sfx):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        hk, bins = histogram(k['k_energy'] / 1e15, bins=logspace(-1, 2,
                                                                 200))
        hh, bins = histogram(h['k_energy'] / 1e15, bins=logspace(-1, 2,
                                                                 200))
        x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
        hr = 1. * hh / hk
        plot(x, hr, label=label)
    xlabel("Energy (PeV)")
    xscale('log')
    ylabel("Ratio")
    legend(loc='best')
    savefig("plots/region_energy_ratio_%s.pdf" % sfx)

def plot_regions_zenith(kevents, hevents, labels, sfx):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(rad2deg(k['k_zenith']), bins=200,
             histtype='step', label=label, normed=True)
        subplot(212)
        hist(rad2deg(h['k_zenith']), bins=200,
             histtype='step', label=label, normed=True)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    legend(loc='best')
    subplot(212)
    xlabel("Zenith $(^\circ)$")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/region_zenith_%s.pdf" % sfx)

    figure()
    for k, h, label in zip(kevents, hevents, labels):
        hk, bins = histogram(rad2deg(k['k_zenith']), bins=200)
        hh, bins = histogram(rad2deg(h['k_zenith']), bins=200)
        x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
        hr = 1. * hh / hk
        plot(x, hr, label=label)
    xlabel("Zenith $(^\circ)$")
    ylabel("Ratio")
    legend(loc='best')
    savefig("plots/region_zenith_ratio_%s.pdf" % sfx)

def plot_regions_azimuth(kevents, hevents, labels, sfx):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(rad2deg(k['k_azimuth']), bins=200,
             histtype='step', label=label, normed=True)
        subplot(212)
        hist(rad2deg(h['k_azimuth']), bins=200,
             histtype='step', label=label, normed=True)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    legend(loc='best')
    subplot(212)
    xlabel("Azimuth $(^\circ)$")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/region_azimuth_%s.pdf" % sfx)

    figure()
    for k, h, label in zip(kevents, hevents, labels):
        hk, bins = histogram(rad2deg(k['k_azimuth']), bins=200)
        hh, bins = histogram(rad2deg(h['k_azimuth']), bins=200)
        x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
        hr = 1. * hh / hk
        plot(x, hr, label=label)
    xlabel("Azimuth $(^\circ)$")
    ylabel("Ratio")
    legend(loc='best')
    savefig("plots/region_azimuth_ratio_%s.pdf" % sfx)

def plot_regions_T200(kevents, hevents, labels, sfx):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(k['k_T200'], bins=200,
             histtype='step', label=label, normed=True)
        subplot(212)
        hist(h['k_T200'], bins=200,
             histtype='step', label=label, normed=True)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    ylim(ymax=.2)
    legend(loc='best')
    subplot(212)
    xlabel("Temperature $(^\circ\mathrm{C})$")
    ylabel("Count")
    ylim(ymax=.2)
    legend(loc='best')
    savefig("plots/region_T200_%s.pdf" % sfx)

    figure()
    for k, h, label in zip(kevents, hevents, labels):
        hk, bins = histogram(k['k_T200'], bins=200)
        hh, bins = histogram(h['k_T200'], bins=200)
        x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
        hr = 1. * hh / hk
        plot(x, hr, label=label)
    xlabel("Temperature $(^\circ\mathrm{C})$")
    ylabel("Ratio")
    legend(loc='best')
    savefig("plots/region_T200_ratio_%s.pdf" % sfx)

def plot_regions_P200(kevents, hevents, labels, sfx):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(k['k_P200'], bins=50,
             histtype='step', label=label, normed=True)
        subplot(212)
        hist(h['k_P200'], bins=50,
             histtype='step', label=label, normed=True)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    ylim(ymax=.2)
    legend(loc='best')
    subplot(212)
    xlabel("Pressure (hPa)")
    ylabel("Count")
    ylim(ymax=.2)
    legend(loc='best')
    savefig("plots/region_P200_%s.pdf" % sfx)

    figure()
    for k, h, label in zip(kevents, hevents, labels):
        hk, bins = histogram(k['k_P200'], bins=50)
        hh, bins = histogram(h['k_P200'], bins=50)
        x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
        hr = 1. * hh / hk
        plot(x, hr, label=label)
    xlabel("Pressure (hPa)")
    ylabel("Ratio")
    legend(loc='best')
    savefig("plots/region_P200_ratio_%s.pdf" % sfx)

def get_regions_data():
    events = data.getNode(GROUP, 'events')

    edges = [0, 15, 30, 300]
    labels = ["%d <= R < %d" % (u, v) for u, v in zip(edges[:-1],
                                                      edges[1:])]
    kevents = []
    hevents = []
    cond = '(r1 <= core_dist) & (core_dist < r2)'
    for r1, r2 in zip(edges[:-1], edges[1:]):
        kevents.append(events.readWhere(cond))
        hevents.append(events.readWhere(cond +
                                        ' & (self_triggered == True)'))
    return kevents, hevents, labels, 'R'

def get_p_regions_data():
    events = data.getNode(GROUP, 'events')

    edges = [995, 1005, 1015]
    labels = ["%d <= P < %d" % (u, v) for u, v in zip(edges[:-1],
                                                      edges[1:])]
    kevents = []
    hevents = []
    cond = '(r1 <= k_P200) & (k_P200 < r2)'
    for r1, r2 in zip(edges[:-1], edges[1:]):
        kevents.append(events.readWhere(cond))
        hevents.append(events.readWhere(cond +
                                        ' & (self_triggered == True)'))
    return kevents, hevents, labels, 'P'

def plot_reconstruction_efficiency_zenith():
    plot_rec_eff('k_theta', bins=linspace(0, deg2rad(40), 50), xf=rad2deg)
    xlabel("KASCADE zenith (deg)")
    savefig('plots/rec_eff_zenith.pdf')

def plot_reconstruction_efficiency_azimuth():
    plot_rec_eff('k_phi', bins=linspace(-pi, pi, 50))
    xlabel("KASCADE azimuth (deg)")
    savefig('plots/rec_eff_azimuth.pdf')
    
    plot_rec_eff('h_phi', bins=linspace(-pi, pi, 50))
    xlabel("HiSPARC azimuth (deg)")
    savefig('plots/rec_eff_azimuth_hisparc.pdf')

def plot_reconstruction_efficiency_density():
    global kevents
    events = data.root.efficiency.events
    #q = '(.5e15 <= k_energy) & (k_energy < 2e15)'
    q = '(k_energy > 2e15)'
    kevents = events.readWhere(q)
    print len(events), len(kevents), 1. * len(kevents) / len(events)
    hevents = events.readWhere(q + '& (self_triggered == True)')
    #hevents = events.readWhere('n_high >= 2')
    #hevents = kevents.compress((kevents['pulseheights'] >= 400.).sum(1) >= 2)
    re1 = kevents.compress((kevents['pulseheights'] >= [20, 0, 20, 20]).all(1))
    re2 = kevents.compress((kevents['pulseheights'] >= [123, 0, 123, 123]).all(1))
    re3 = kevents.compress((kevents['pulseheights'] >= [210, 0, 210, 210]).all(1))

    rer1 = re1.compress(re1['reconstructed'] == True)
    rer2 = re2.compress(re2['reconstructed'] == True)

    figure()
    bins = linspace(0, 10, 50)
    DD = .37
    plot_hist_ratio(hevents['k_cosdens_e'] + DD, kevents['k_cosdens_e'] + DD,
                    bins=bins, label="trigger")
    plot_hist_ratio(re1['k_cosdens_e'] + DD, kevents['k_cosdens_e'] + DD,
                    bins=bins, label="20 ADC")
    #plot_hist_ratio(rer1['k_cosdens_e'], kevents['k_cosdens_e'],
    #                bins=linspace(0, 10, 50), label="20 ADC reconstructed")
    plot_hist_ratio(re2['k_cosdens_e'] + DD, kevents['k_cosdens_e'] + DD,
                    bins=bins, label="70 mV in corners")
    plot_hist_ratio(re3['k_cosdens_e'] + DD, kevents['k_cosdens_e'] + DD,
                    bins=bins, label="120 mV in corners")
    #plot_hist_ratio(rer2['k_cosdens_e'] + DD, kevents['k_cosdens_e'] + DD,
    #                bins=linspace(0, 10, 50), label="70 mV reconstructed")

    x = bins
    p0 = exp(-.5 * x)           # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    ptrig = 1 - (pd0 + 4 * pd1)
    pcorn = pd3 + pd4
    #pcorn = 1 - (pd0 + 4 * pd1 + 6 * pd2 + 3 * pd3)  # for testing
    plot(x, ptrig, label="Poisson trigger")
    plot(x, pcorn, label="Poisson corners")

    xlabel("verschoven KASCADE density")
    ylabel("Probability")
    legend(loc='best')
    ylim(0, 1.05)
    title(q)
    savefig('plots/rec_eff_density-trigger.pdf')

    figure()
    plot_hist_ratio(rer1['k_dens_e'], re1['k_dens_e'],
                    bins=linspace(0, 10, 50), label="20 ADC in corners")
    plot_hist_ratio(rer2['k_dens_e'], re2['k_dens_e'],
                    bins=linspace(0, 10, 50), label="70 mV in corners")
    x = linspace(1, 10, 50)
    plot(x, std_t(x), label="Timing uncertainty")
    xlabel("KASCADE density")
    ylabel("Reconstruction efficiency")
    legend(loc='best')
    ylim(0, 1.05)
    savefig('plots/rec_eff_density-reconstruction.pdf')

def plot_reconstruction_efficiency_coredist():
    plot_rec_eff('core_dist', bins=linspace(0, 200, 200), dens=1.)
    xlabel("Core distance [m]")
    savefig('plots/rec_eff_coredist.pdf')

def plot_rec_eff(var, bins, dens=1., xf=None):
    events = data.root.efficiency.events

    kevents = events[:]
    hevents = events.readWhere('self_triggered == True')
    revents = events.readWhere('reconstructed == True')

    kevents = kevents.compress(kevents['k_dens_e'][:,1] >= dens)
    hevents = hevents.compress(hevents['k_dens_e'][:,1] >= dens)
    revents = revents.compress(revents['k_dens_e'][:,1] >= dens)

    figure()
    plot_hist_ratio(hevents[var], kevents[var], bins, xf=xf,
                    label="trigger efficiency")
    plot_hist_ratio(revents[var], hevents[var], bins, xf=xf,
                    label="reconstruction efficiency")
    plot_hist_ratio(revents[var], kevents[var], bins, xf=xf,
                    label="combined efficiency")
    ylim(0, 1.05)
    ylabel("Count")
    title("Electron density >= %.1f" % dens)
    legend(loc='best')

def plot_hist_ratio(events1, events2, bins, xf=None, label=None):
    x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
    if xf:
        x = xf(x)
    h1, bins = histogram(events1, bins=bins)
    h2, bins = histogram(events2, bins=bins)
    hr = 1. * h1 / h2
    plot(x, hr, label=label)

def plot_detector_efficiency():
    pass 

def plot_rec_eff_timings():
    global sel, dt1, dt2, m, r, rmin
    events = data.root.efficiency.events


    clf()

    sel = events.readWhere('reconstructed == True')
    r = array([rvalue_from_event(u) for u in progressbar(sel)])
    r = r.compress(-isnan(r))
    hist(r, bins=arange(0, 2, .1), histtype='step', label="reconstructed")

    sel = events.readWhere('reconstructed == False')
    r = array([rvalue_from_event(u) for u in progressbar(sel)])
    r = r.compress(-isnan(r))
    hist(r, bins=arange(0, 10, .1), histtype='step',
         label="max, not reconstructed")

    xlabel("Fraction of station size traversed")
    ylabel("Count")
    legend()
    savefig("plots/rec_eff_timings.pdf")

def rvalue_from_event(event):
    t1, t2, t3, t4 = event['t']
    ph1, ph2, ph3, ph4 = event['pulseheights']

    dt1 = t1 - t3
    dt2 = t1 - t4
    if isnan(dt1) or isnan(dt2):
        return nan
    elif min(ph1, ph3, ph4) <= 70 / .57:
        return nan
    else:
        x1 = 10 * cos(event['h_phi'] - phi1)
        x2 = 10 * cos(event['h_phi'] - phi2)
        r1 = 3e-1 * dt1 / x1
        r2 = 3e-1 * dt2 / x2
        if r1 != 0.:
            return r1
        else:
            return r2
    raise RuntimeError("Should not reach this point")

def plot_rvalue_core():
    events = data.root.efficiency.events
    q = "(.8e15 <= k_energy) & (k_energy < 2e15)"

    global x, y, yerr, y2
   
    figure()

    bins = linspace(0, 100, 21)
    x, y, yerr = [], [], []
    y2 = []
    #for r0, r1 in progressbar(zip(bins[:-1], bins[1:])):
    print "(r0\tr1)\tcount\t>20 ADC\t>30 mV\t>70 mV\t>120 mV"
    for r0, r1 in (zip(bins[:-1], bins[1:])):
        sel = events.readWhere(q + '& (r0 <= core_dist) & (core_dist < r1)')
        c, c20adc, c30mV, c70mV, c120mV = 0, 0, 0, 0, 0
        for event in sel:
            ph1, ph2, ph3, ph4 = event['pulseheights']
            level = min(ph1, ph3, ph4)
            c += 1
            if level >= 20:
                c20adc += 1
            if level >= 30 / .57:
                c30mV += 1
            if level >= 70 / .57:
                c70mV += 1
            if level >= 120 / .57:
                c120mV += 1
        print "(%.1f\t%.1f)\t%d\t%d\t%d\t(%.2f)\t%d\t(%.2f)\t%d\t(%.2f)" % (r0, r1, c, c20adc, c30mV, c30mV / c20adc, c70mV, c70mV / c20adc, c120mV, c120mV / c20adc)

        #x.append(mean([r0, r1]))
        #rvs = array([rvalue_from_event(u) for u in sel])
        #rvs = rvs.compress(-isnan(rvs))
        #y.append(median(rvs))
        #y2.append(rvs)
        #yerr.append(std(rvs))
        #print 1. * len(rvs) / len(sel), len(rvs)

    #errorbar(x, y, yerr=yerr)
    plot(x, y)
    xlabel("Core distance [m]")
    ylabel("rvalue")
    savefig("plots/rvalues-core-dist.pdf")

def electron_gamma_trigger():
    events = data.root.efficiency.events
    q = '(.5e15 <= k_energy) & (k_energy < 2e15)'
    kevents = events.readWhere(q)
    hevents = events.readWhere(q + '& (self_triggered == True)')
    re1 = kevents.compress((kevents['pulseheights'] >= [300, 300, 300, 300]).sum(1) >= 2)

    clf()
    bins = linspace(0, 10, 50)
    DD = 0
    plot_hist_ratio(hevents['k_cosdens_e'] + DD, kevents['k_cosdens_e'] + DD,
                    bins=bins, label="trigger")
    plot_hist_ratio(re1['k_cosdens_e'] + re1['k_dens_mu'] + DD,
                    kevents['k_cosdens_e'] + kevents['k_dens_mu'] + DD,
                    bins=bins, label="300 ADC")

    x = bins
    p0 = exp(-.5 * x)           # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    ptrig = 1 - (pd0 + 4 * pd1)
    plot(x, ptrig, label="Poisson trigger")

    xlabel("KASCADE density")
    ylabel("Probability")
    legend(loc='best')
    ylim(0, 1.05)
    title(q)
    savefig('plots/electron_gamma_trigger.pdf')

def plot_poisson_trigger_cuts():
    events = data.getNode(GROUP, 'events')

    #q = '(k_Num_e > 10 ** 4.6) & (k_zenith <= %f) & (core_dist_center <= 90)' % deg2rad(30)
    q = '(k_Num_e > 10 ** 4.6) & (k_Num_e < 10 ** 7) & (k_zenith <= %f) & (core_dist_center <= 90)' % deg2rad(30)
    #q = '(k_Num_e > 0)'

    clf()
    bins = linspace(0, 10, 201)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    kevents = events.readWhere(q)
    hevents = events.readWhere(q + '& (self_triggered == True)')
    kdens = mean(kevents['k_dens_e'] + kevents['k_dens_mu'], 1) * \
            cos(kevents['k_theta'])
    hdens = mean(hevents['k_dens_e'] + hevents['k_dens_mu'], 1) * \
            cos(hevents['k_theta'])

    hk, bins = histogram(kdens, bins=bins)
    hh, bins = histogram(hdens, bins=bins)
    hr = 1. * hh / hk
    plot(x, hr, label="Data")
    
    # Poisson probabilities
    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    ptrig = 1 - (pd0 + 4 * pd1)         # Trigger probability
    plot(x, ptrig, label="Poisson")

    xlabel("Electron density (m$^{-1}$)")
    ylabel("Probability of trigger")
    legend(loc='best')
    xlim(0, 10)
    ylim(0, 1.05)
    savefig("plots/poisson_trigger_cuts.pdf")

def plot_charged_particles_poisson(use_known=False):
    events = data.getNode(GROUP, 'events')
    s = Scintillator()

    if use_known:
        s.mev_scale = 0.0086306040898338834
        s.gauss_scale = 0.84289265239940525
    else:
        ph = events[:]['pulseheights'][:,1]
        analyze_charged_particle_spectrum(s, ph, constrained=False)

    bins = linspace(0, 5, 41)
    x = bins[:-1] + .5 * (bins[1] - bins[0])
    y, yerr = [], []

    events = events.read()
    dens = events['k_cosdens_charged'][:,1]     # center detector
    for num, (u, v) in enumerate(zip(bins[:-1], bins[1:])):
        print "Analyzing %.2f <= dens < %.2f" % (u, v)
        sel = events.compress((u <= dens) & (dens < v))
        ph = sel['pulseheights'][:,1]
        p = analyze_charged_particle_spectrum(s, ph, constrained=True)
        title("Charged particle part of spectrum (%.2f <= dens < %.2f)" %
              (u, v))
        savefig("plots/charged_particles_poisson_bin_%d.pdf" % num)
        y.append(p[0])
        yerr.append(p[1])
    y = array(y)
    yerr = array(yerr)

    figure()
    errorbar(x, y, yerr=yerr, fmt='o', label="Data")

    # Poisson probabilities
    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # one or more particles

    plot(x, pp, label="Poisson")

    xlabel("Charged particle density (m$^{-2}$)")
    ylabel("Probability of one or more particles")
    legend(loc='best')
    savefig("plots/charged_particles_poisson.pdf")

    if 'poisson' in data.root.datasets:
        data.removeNode('/datasets', 'poisson', recursive=True)

    data.createGroup('/datasets', 'poisson')
    data.createArray('/datasets/poisson', 'x', x)
    data.createArray('/datasets/poisson', 'y', y)
    data.createArray('/datasets/poisson', 'yerr', yerr)

def analyze_charged_particle_spectrum(s, ph, constrained=False):
    if not constrained:
        x0 = (10 ** 4, 3.38 / 380., 1)
        residuals = s.residuals
    else:
        x0 = 10
        residuals = s.constrained_residuals

    figure()
    # Fit of convoluted Landau
    n, bins = histogram(ph, bins=linspace(0, 2000, 101))
    nx = bins[:-1] + .5 * (bins[1] - bins[0])
    x = linspace(-2000, 2000, 201)
    y = interp(x, nx, n)
    xopt, fopt, iter, funcalls, warnflag = fmin(residuals, x0,
                                                (x, y, 350, 500), disp=0,
                                                full_output=1)
    plot(x, s.conv_landau(x, *xopt))

    print "Residuals: %.2f" % fopt

    # Charged particle spectrum
    step(x, y, where='mid')
    yl = s.conv_landau(x, *xopt)
    plot(x, yl)
    i = (y <= yl).argmax()
    yp = array(yl[:i].tolist() + y[i:].tolist())
    step(x, yp, where='mid')
    N_T = sum(y.compress(x >= 0))
    N_CP = sum(yp.compress(x >= 0))
    print "Charged particles: %.2f %% of events" % ((N_CP / N_T) * 100)

    xlim(xmin=0)
    yscale('log')
    ylim(ymin=1)
    xlabel("Pulseheight [ADC counts]")
    ylabel("Counts")

    return N_CP / N_T, sqrt(N_CP) / N_T

def plot_conv_poisson():
    figure()

    x = data.root.datasets.poisson.x.read()
    y = data.root.datasets.poisson.y.read()
    yerr = data.root.datasets.poisson.yerr.read()
    errorbar(x, y, yerr=yerr, fmt='o', label="Data")

    x = linspace(-10, 10, 101)
    f = vectorize(lambda x: 1 - exp(-.5 * x) if x >= 0 else 0.)

    for s in [0., 0.5, 1., 1.5]:
        if s == 0.:
            plot(x, f(x), label="Unconvoluted")
        else:
            g = stats.norm(scale=s).pdf
            plot(x, discrete_convolution(f, g, x),
                 label="sigma = %.2f m$^{-2}$" % s)

    xlim(0, 5)
    title("Effect of uncertainty on particle density")
    xlabel("Charged particle density (m$^{-2}$)")
    ylabel("Probability of one or more particles")
    legend(loc='best')
    savefig("plots/conv_poisson.pdf")


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'a')

    main()
