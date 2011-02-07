import tables

from build_dataset import DETECTORS
from plot_utils import *

from pylab import *
from scipy.optimize import curve_fit


DATAFILE = 'kascade.h5'
GROUP = '/efficiency'


def main():
    #plot_energy_histogram()
    #plot_core_distance()
    #plot_core_alpha()
    #plot_core_pos()
    #plot_zenith()
    #plot_azimuth()
    #plot_T200()
    #plot_P200()
    #plot_particle_density()
    plot_density_fraction_err()
    #plot_density_fraction()
    #optimize_trigger_prob()
    #plot_theta_phi_corr()

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
    events = data.getNode(GROUP, 'events')

    global err, hdens

    figure()
    bins = linspace(0, 10, 201)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    kevents = events[:]
    hevents = events.readWhere('self_triggered == True')
    kdens = kevents['k_dens_e']
    kdens_m = median(kdens, 1)
    hdens = hevents['k_dens_e']
    hdens_m = median(hdens, 1)

    err = [[] for u in xrange(len(bins) - 1)]
    for idx, v in zip(bins.searchsorted(hdens_m), hdens):
        try:
            err[idx - 1].extend(v.tolist())
        except IndexError:
            pass
    err = [std(u) for u in err]

    hk, bins = histogram(kdens_m, bins=bins)
    hh, bins = histogram(hdens_m, bins=bins)
    hr = 1. * hh / hk
    errorbar(x, hr, xerr=err, label="Data", capsize=0)
    
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
    events = data.getNode(GROUP, 'events')

    figure()
    bins = linspace(0, 10, 201)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    for Z in ['k_zenith < %.2f' % deg2rad(40), 'k_zenith < %.2f' %
              deg2rad(10), 'k_zenith > %.2f' % deg2rad(25)]:
        kevents = events.readWhere(Z)
        hevents = events.readWhere('(%s) & (self_triggered == True)' % Z)
        kdens = median(kevents['k_dens_e'], 1)
        hdens = median(hevents['k_dens_e'], 1)

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

    global hrs, x, p

    figure()
    bins = linspace(0, 10, 201)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
    hrs = []
    for Z in ['k_zenith < %.2f' % deg2rad(40), 'k_zenith < %.2f' %
              deg2rad(10), 'k_zenith > %.2f' % deg2rad(25)]:
        kevents = events.readWhere(Z)
        hevents = events.readWhere('(%s) & (self_triggered == True)' % Z)
        kdens = median(kevents['k_dens_e'], 1) * cos(kevents['k_zenith'])
        hdens = median(hevents['k_dens_e'], 1) * cos(hevents['k_zenith'])

        hk, bins = histogram(kdens, bins=bins)
        hh, bins = histogram(hdens, bins=bins)
        hr = 1. * hh / hk
        hrs.append(hr)

    # Poisson probability of zero particles in detector
    p = lambda x, c: 1 - (4 * (1 - exp(-.5 * (x + c))) * \
                     exp(-.5 * (x + c)) ** 3 + exp(-.5 * (x + c)) ** 4)

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
    global events
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



if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    main()
