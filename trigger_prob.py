import tables

from pylab import *


#DATAFILE = 'trigprob.h5'
DATAFILE = 'data-e15.h5'
KASCADEDATAFILE = 'kascade.h5'


def main():
    plot_dens_from_station_events()
    plot_trigger_probability()
    plot_density_spread()
    plot_density_coredist()

def plot_dens_from_station_events():
    events = data.root.analysis.angle_0

    global n

    x, y = [], []
    for event in events:
        n = array([event['n1'], event['n2'], event['n3'], event['n4']],
                  dtype=float64)
        n += uniform(-.5, .5, 4)
        n *= 2.

        x.append(mean(n))
        y.append(max(n) - min(n))

    figure()
    plot(x, y, ',')
    xlabel("Mean density")
    ylabel("Density spread (max - min)")
    savefig('plots/trigprob-dens-station-events.pdf')

def plot_trigger_probability():
    global events, dens, bins, x, trig, all

    events = data.root.analysis.angle_0
    k_events = kdata.getNode('/efficiency', 'events')
    q = '(.9e15 <= k_energy) & (k_energy < 1.1e15) & (k_zenith < .35)'
    kevents = k_events.readWhere(q)
    ktevents = k_events.readWhere(q + '& (self_triggered == True)')
   
    bins = linspace(0, 10., 101)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']
    kdens = array([mean(u['k_dens_e']) for u in kevents[:]])
    ktdens = array([mean(u['k_dens_e']) for u in ktevents[:]])

    all = zeros(len(bins) - 1)
    trig = zeros(len(bins) - 1)
    for idx, event in zip(bins.searchsorted(dens), events[:]):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        if sum([int(event[u] >= 1.) for u in ['n1', 'n2', 'n3', 'n4']]) >= 2:
            trig[idx - 1] += 1

    kh, bins = histogram(kdens, bins=bins)
    kth, bins = histogram(ktdens, bins=bins)
    trigprob = 1. * kth / kh

    figure()
    y = trig / all
    plot(x, y, label="Simulation")
    plot(x, trigprob, label="KASCADE")

    # Poisson probability of zero particles in detector
    x = linspace(0, 10, 101)
    p0 = exp(-.5 * x)
    p = 1 - (4 * (1 - p0) * p0 ** 3 + p0 ** 4)
    plot(x, p, label="Poisson")

    xlabel("Mean density")
    ylabel("Trigger probability")
    legend(loc='best')
    ylim(0, 1.05)
    title(q)
    savefig('plots/trigprob-trigger-probability.pdf')

def plot_density_spread():
    events = data.root.analysis.angle_0
    kevents = kdata.getNode('/efficiency', 'events')
    q = '(.9e15 <= k_energy) & (k_energy < 1.1e15) & (k_zenith < .35)'
    kevents = kevents.readWhere(q)

    dmean = array([mean([u['n1'], u['n2'], u['n3'], u['n4']]) * 2. for u
                   in events[:]])
    dspread = array([(max([u['n1'], u['n2'], u['n3'], u['n4']]) -
                     min([u['n1'], u['n2'], u['n3'], u['n4']])) * 2. for u
                     in events[:]])
    kdmean = array([mean(u['k_dens_e']) for u in kevents[:]])
    kdspread = array([max(u['k_dens_e']) - min(u['k_dens_e']) for u in
                     kevents[:]])

    bins = arange(0, 10.1, .5)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    xspread = [compress((u <= dmean) & (dmean < v), dspread) for
                     u, v in zip(bins[:-1], bins[1:])]
    y = array([mean(u) for u in xspread])
    yerr = array([std(u) for u in xspread])

    kxspread = [compress((u <= kdmean) & (kdmean < v), kdspread) for
                     u, v in zip(bins[:-1], bins[1:])]
    ky = array([mean(u) for u in kxspread])
    kyerr = array([std(u) for u in kxspread])

    figure()
    errorbar(x, y, yerr=yerr, capsize=0, label="Simulation")
    errorbar(x, ky, yerr=kyerr, capsize=0, label="KASCADE")
    legend(loc='best')
    xlabel("Mean density (m$^{-1}$)")
    ylabel("Density spread (m$^{-1}$)")
    title(q)
    savefig('plots/trigprob-density-spread.pdf')

def plot_density_coredist():
    events = data.root.analysis.angle_0
    kevents = kdata.getNode('/efficiency', 'events')
    q = '(.9e15 <= k_energy) & (k_energy < 1.1e15) & (k_zenith < .35)'
    kevents = kevents.readWhere(q)

    dmean = array([mean([u['n1'], u['n2'], u['n3'], u['n4']]) * 2. for u
                   in events[:]])
    coredist = events[:]['r']

    kdmean = array([mean(u['k_dens_e']) for u in kevents[:]])
    kcoredist = kevents[:]['core_dist']

    bins = arange(0, 10.1, .5)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    xcoredist = [compress((u <= dmean) & (dmean < v), coredist) for u, v
                 in zip(bins[:-1], bins[1:])]
    y = array([mean(u) for u in xcoredist])
    yerr = array([std(u) for u in xcoredist])

    kxcoredist = [compress((u <= kdmean) & (kdmean < v), kcoredist) for u, v
                 in zip(bins[:-1], bins[1:])]
    ky = array([mean(u) for u in kxcoredist])
    kyerr = array([std(u) for u in kxcoredist])


    figure()
    errorbar(x, y, yerr=yerr, capsize=0, label="Simulation")
    errorbar(x, ky, yerr=kyerr, capsize=0, label="KASCADE")
    legend(loc='best')
    xlabel("Mean density (m$^{-1}$)")
    ylabel("Mean core distance (m)")
    title(q)
    savefig('plots/trigprob-density-coredist.pdf')


if __name__ == '__main__':
    try:
        data, kdata
    except NameError:
        data = tables.openFile(DATAFILE, 'r')
        kdata = tables.openFile(KASCADEDATAFILE, 'r')

    main()
