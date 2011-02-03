import tables

from build_dataset import DETECTORS
from plot_utils import *

from pylab import *


DATAFILE = 'kascade.h5'
GROUP = '/efficiency'


def main():
    plot_energy_histogram()
    plot_core_distance()
    plot_core_alpha()
    plot_core_pos()
    plot_zenith()
    plot_azimuth()
    plot_T200()
    plot_P200()

    kevents, hevents, labels = get_regions_data()
    plot_regions_core_pos(kevents, hevents, labels)
    plot_regions_energy_count(kevents, hevents, labels)
    plot_regions_energy_ratio(kevents, hevents, labels)
    plot_regions_zenith(kevents, hevents, labels)
    plot_regions_azimuth(kevents, hevents, labels)
    plot_regions_T200(kevents, hevents, labels)
    plot_regions_P200(kevents, hevents, labels)

def plot_energy_histogram():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['k_energy'] / 1e15, bins=logspace(-1, 2, 200),
         histtype='step', label="KASCADE trigger")
    hist(events.readWhere('self_triggered == True')['k_energy'] / 1e15,
         bins=logspace(-1, 2, 200), histtype='step',
         label="HiSPARC trigger")
    xlabel("Energy (PeV)")
    xscale('log')
    ylabel("Count")
    legend(loc='best')
    savefig("plots/energy_histogram.pdf")

def plot_core_distance():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['core_dist'], bins=linspace(0, 200, 200),
         histtype='step', label="KASCADE trigger")
    hist(events.readWhere('self_triggered == True')['core_dist'],
         bins=linspace(0, 200, 200), histtype='step',
         label="HiSPARC trigger")
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
         label="KASCADE trigger")
    hist(events.readWhere('self_triggered == True')['core_alpha'] / pi,
         bins=200, histtype='step', label="HiSPARC trigger")
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
         label="KASCADE trigger")
    hist(rad2deg(events.readWhere('self_triggered == True')['k_zenith']),
         bins=200, histtype='step', label="HiSPARC trigger")
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
         label="KASCADE trigger")
    hist(rad2deg(events.readWhere('self_triggered == True')['k_azimuth']),
         bins=200, histtype='step', label="HiSPARC trigger")
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
    hist(events[:]['k_T200'], bins=200, histtype='step',
         label="KASCADE trigger")
    hist(events.readWhere('self_triggered == True')['k_T200'], bins=200,
         histtype='step', label="HiSPARC trigger")
    xlabel("Temperature $(^\circ\!\mathrm{C})$")
    ylabel("Count")
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
         label="KASCADE trigger")
    hist(events.readWhere('self_triggered == True')['k_P200'], bins=200,
         histtype='step', label="HiSPARC trigger")
    xlabel("Pressure (hPa)")
    ylabel("Count")
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

def plot_regions_core_pos(kevents, hevents, labels):
    for i, (k, h, label) in enumerate(zip(kevents, hevents, labels)):
        figure()
        cores = k[:]['k_core_pos']
        contour_histogram2d(cores[:,0], cores[:,1], bins=200)
        plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
        xlabel("(m)")
        ylabel("(m)")
        title("Shower core positions (%s)" % label)
        savefig("plots/region_core_pos_%d.pdf" % i)

    for i, (k, h, label) in enumerate(zip(kevents, hevents, labels)):
        figure()
        cores = h[:]['k_core_pos']
        contour_histogram2d(cores[:,0], cores[:,1], bins=200)
        plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
        xlabel("(m)")
        ylabel("(m)")
        title("Shower core positions (self-triggered) (%s)" % label)
        savefig("plots/region_core_pos_trig_%d.pdf" % i)

def plot_regions_energy_count(kevents, hevents, labels):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(k['k_energy'] / 1e15, bins=logspace(-1, 2, 200),
             histtype='step', label=label)
        subplot(212)
        hist(h['k_energy'] / 1e15, bins=logspace(-1, 2, 200),
             histtype='step', label=label)
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
    savefig("plots/region_energy.pdf")

def plot_regions_energy_ratio(kevents, hevents, labels):
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
    savefig("plots/region_energy_ratio.pdf")

def plot_regions_zenith(kevents, hevents, labels):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(rad2deg(k['k_zenith']), bins=200,
             histtype='step', label=label)
        subplot(212)
        hist(rad2deg(h['k_zenith']), bins=200,
             histtype='step', label=label)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    legend(loc='best')
    subplot(212)
    xlabel("Zenith $(^\circ)$")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/region_zenith.pdf")

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
    savefig("plots/region_zenith_ratio.pdf")

def plot_regions_azimuth(kevents, hevents, labels):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(rad2deg(k['k_azimuth']), bins=200,
             histtype='step', label=label)
        subplot(212)
        hist(rad2deg(h['k_azimuth']), bins=200,
             histtype='step', label=label)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    legend(loc='best')
    subplot(212)
    xlabel("Azimuth $(^\circ)$")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/region_azimuth.pdf")

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
    savefig("plots/region_azimuth_ratio.pdf")

def plot_regions_T200(kevents, hevents, labels):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(k['k_T200'], bins=200,
             histtype='step', label=label)
        subplot(212)
        hist(h['k_T200'], bins=200,
             histtype='step', label=label)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    legend(loc='best')
    subplot(212)
    xlabel("Temperature $(^\circ\mathrm{C})$")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/region_T200.pdf")

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
    savefig("plots/region_T200_ratio.pdf")

def plot_regions_P200(kevents, hevents, labels):
    figure()
    for k, h, label in zip(kevents, hevents, labels):
        subplot(211)
        hist(k['k_P200'], bins=200,
             histtype='step', label=label)
        subplot(212)
        hist(h['k_P200'], bins=200,
             histtype='step', label=label)
    subplot(211)
    title("KASCADE (upper), self-triggered (lower)")
    ylabel("Count")
    legend(loc='best')
    subplot(212)
    xlabel("Pressure (hPa)")
    ylabel("Count")
    legend(loc='best')
    savefig("plots/region_P200.pdf")

    figure()
    for k, h, label in zip(kevents, hevents, labels):
        hk, bins = histogram(k['k_P200'], bins=200)
        hh, bins = histogram(h['k_P200'], bins=200)
        x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
        hr = 1. * hh / hk
        plot(x, hr, label=label)
    xlabel("Pressure (hPa)")
    ylabel("Ratio")
    legend(loc='best')
    savefig("plots/region_P200_ratio.pdf")

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
    return kevents, hevents, labels



if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    main()
