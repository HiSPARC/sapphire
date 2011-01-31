import tables

from build_dataset import DETECTORS

from pylab import *


DATAFILE = 'kascade.h5'
GROUP = '/efficiency'


def main():
    #plot_energy_histogram()
    plot_core_distance()
    #plot_core_alpha()
    #plot_core_pos()

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
    legend()
    savefig("plots/energy_histogram.pdf")

def plot_core_distance():
    events = data.getNode(GROUP, 'events')

    figure()
    hist(events[:]['core_dist'], bins=200, histtype='step',
         label="KASCADE trigger")
    hist(events.readWhere('self_triggered == True')['core_dist'],
         bins=200, histtype='step', label="HiSPARC trigger")
    xlabel("Core distance (m)")
    ylabel("Count")
    legend()
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
    legend()
    savefig("plots/core_alpha.pdf")

def plot_core_pos():
    events = data.getNode(GROUP, 'events')

    mylog = vectorize(lambda x: log10(x) if x > 0 else 0.)

    figure()
    cores = events[:]['k_core_pos']
    H, xedges, yedges = histogram2d(cores[:,0], cores[:,1], bins=200)
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2
    contourf(x, y, H.T, 10)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions")
    savefig("plots/core_pos.pdf")

    figure()
    H = mylog(H)
    contourf(x, y, H.T, 10)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions")
    savefig("plots/core_pos_log.pdf")
    
    figure()
    cores = events.readWhere('self_triggered == True')['k_core_pos']
    H, xedges, yedges = histogram2d(cores[:,0], cores[:,1], bins=200)
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2
    contourf(x, y, H.T, 10)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions")
    savefig("plots/core_pos_trig.pdf")

    figure()
    H = mylog(H)
    contourf(x, y, H.T, 10)
    plot([u[0] for u in DETECTORS], [u[1] for u in DETECTORS], 'wo')
    xlabel("(m)")
    ylabel("(m)")
    title("Shower core positions")
    savefig("plots/core_pos_trig_log.pdf")


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    main()
