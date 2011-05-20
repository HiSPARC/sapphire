import tables
import progressbar as pb

from pylab import *


#DATAFILE = 'trigprob.h5'
DATAFILE = 'data-e15.h5'
KASCADEDATAFILE = 'kascade.h5'


def progressbar(*args, **kwargs):
    progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])
    return progress(*args, **kwargs)

def main():
    #plot_dens_from_station_events()
    #plot_trigger_probability()
    plot_reconstruction_probability()
    #plot_poisson_probability()
    #plot_density_spread()
    #plot_density_coredist()
    #plot_poisson_zero()
    #plot_poisson_one()
    #plot_poisson_two()
    #plot_poisson_three()
    #plot_poisson_four()

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
    for idx, event in progressbar(zip(bins.searchsorted(dens), events[:])):
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

def plot_reconstruction_probability():
    global d34, all

    events = data.root.analysis.angle_0
   
    bins = linspace(0, 10, 51)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']

    all = zeros(len(bins) - 1)
    d34 = zeros(len(bins) - 1)
    for idx, event in progressbar(zip(bins.searchsorted(dens), events[:])):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        if event['n1'] >= 1 and event['n2'] == 0 and event['n3'] >= 1 \
            and event['n4'] >= 1:
            d34[idx - 1] += 1
        elif event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] >= 1 \
            and event['n4'] >= 1:
            d34[idx - 1] += 1

    clf()

    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    plot(x, pd3 + pd4, label="Poisson")
    plot(x, d34 / all, label="Simulation")

    xlabel("Mean particle density [m$^{-2}$]")
    ylabel("Probability")
    legend(loc='best', numpoints=1)
    ylim(0, 1.05)
    title("Probabilities of number of stations hit")
    savefig('plots/trigprob-reconstruction-probability.pdf')

def plot_poisson_probability():
    global events, dens, bins, x, d1, all

    events = data.root.analysis.angle_0
   
    bins = linspace(0, 10, 51)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']

    all = zeros(len(bins) - 1)
    d0 = zeros(len(bins) - 1)
    d1 = zeros(len(bins) - 1)
    d2 = zeros(len(bins) - 1)
    d3 = zeros(len(bins) - 1)
    d4 = zeros(len(bins) - 1)
    for idx, event in progressbar(zip(bins.searchsorted(dens), events[:])):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        n = sum([event[u] >= 1 for u in ['n1', 'n2', 'n3', 'n4']])
        if n == 0:
            d0[idx - 1] += 1
        elif n == 1:
            d1[idx - 1] += 1
        elif n == 2:
            d2[idx - 1] += 1
        elif n == 3:
            d3[idx - 1] += 1
        elif n == 4:
            d4[idx - 1] += 1

        #if event['n1'] == 0 and event['n2'] == 0 and event['n3'] == 0 \
        #    and event['n4'] == 0:
        #    d0[idx - 1] += 1
        #if event['n1'] >= 1 and event['n2'] == 0 and event['n3'] == 0 \
        #    and event['n4'] == 0:
        #    d1[idx - 1] += 1
        #if event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] == 0 \
        #    and event['n4'] == 0:
        #    d2[idx - 1] += 1
        #if event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] >= 1 \
        #    and event['n4'] == 0:
        #    d3[idx - 1] += 1
        #if event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] >= 1 \
        #    and event['n4'] >= 1:
        #    d4[idx - 1] += 1

    clf()

    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    plot(x, pd0, label="calc: 0 hit")
    plot(x, d0 / all, 'o', label="sim: 0 hit")
    plot(x, pd1, label="calc: 1 hit")
    plot(x, d1 / 4 / all, 'o', label="sim: 1 hit")
    plot(x, pd2, label="calc: 2 hit")
    plot(x, d2 / 6 / all, 'o', label="sim: 2 hit")
    plot(x, pd3, label="calc: 3 hit")
    plot(x, d3 / 4 / all, 'o', label="sim: 3 hit")
    plot(x, pd4, label="calc: 4 hit")
    plot(x, d4 / all, 'o', label="sim: 4 hit")

    xlabel("Mean particle density [m$^{-2}$]")
    ylabel("Probability")
    legend(loc='best', numpoints=1)
    ylim(0, 1.05)
    title("Probabilities of number of stations hit")
    savefig('plots/trigprob-poisson-probability.pdf')

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

def plot_poisson_zero():
    global all, d0

    events = data.root.analysis.angle_0
   
    bins = linspace(0, 10, 51)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']

    all = zeros(len(bins) - 1)
    d0 = zeros(len(bins) - 1)
    for idx, event in progressbar(zip(bins.searchsorted(dens),
                                      events[:])):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        if event['n1'] == 0 and event['n2'] == 0 and event['n3'] == 0 \
            and event['n4'] == 0:
            d0[idx - 1] += 1

    clf()

    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    plot(x, pd0, label="Poisson")
    plot(x, d0 / all, label="Data")

    xlabel("Mean density")
    ylabel("Probability")
    title("Prob: no detectors with particles")
    legend(loc='best')
    savefig("plots/poisson-zero.pdf")

def plot_poisson_one():
    events = data.root.analysis.angle_0
   
    bins = linspace(0, 10, 51)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']

    all = zeros(len(bins) - 1)
    d1 = zeros(len(bins) - 1)
    d2 = zeros(len(bins) - 1)
    d3 = zeros(len(bins) - 1)
    d4 = zeros(len(bins) - 1)
    for idx, event in progressbar(zip(bins.searchsorted(dens),
                                      events[:])):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        if event['n1'] >= 1 and event['n2'] == 0 and event['n3'] == 0 \
            and event['n4'] == 0:
            d1[idx - 1] += 1
        if event['n1'] == 0 and event['n2'] >= 1 and event['n3'] == 0 \
            and event['n4'] == 0:
            d2[idx - 1] += 1
        if event['n1'] == 0 and event['n2'] == 0 and event['n3'] >= 1 \
            and event['n4'] == 0:
            d3[idx - 1] += 1
        if event['n1'] == 0 and event['n2'] == 0 and event['n3'] == 0 \
            and event['n4'] >= 1:
            d4[idx - 1] += 1

    figure()

    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    plot(x, pd1, label="Poisson")
    #plot(x, d1 / all, label="Data: 1")
    #plot(x, d2 / all, label="Data: 2")
    #plot(x, d3 / all, label="Data: 3")
    #plot(x, d4 / all, label="Data: 4")
    plot(x, (d1 + d2 + d3 + d4) / 4 / all, label="Simulation")

    xlabel("Mean density")
    ylabel("Probability")
    title("Prob: one detector with particles")
    legend(loc='best')
    savefig("plots/poisson-one.pdf")

def plot_poisson_two():
    events = data.root.analysis.angle_0
   
    bins = linspace(0, 10, 51)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']

    all = zeros(len(bins) - 1)
    d12 = zeros(len(bins) - 1)
    d13 = zeros(len(bins) - 1)
    d14 = zeros(len(bins) - 1)
    d23 = zeros(len(bins) - 1)
    d24 = zeros(len(bins) - 1)
    d34 = zeros(len(bins) - 1)
    for idx, event in progressbar(zip(bins.searchsorted(dens),
                                      events[:])):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        if event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] == 0 \
            and event['n4'] == 0:
            d12[idx - 1] += 1
        if event['n1'] >= 1 and event['n2'] == 0 and event['n3'] >= 1 \
            and event['n4'] == 0:
            d13[idx - 1] += 1
        if event['n1'] >= 1 and event['n2'] == 0 and event['n3'] == 0 \
            and event['n4'] >= 1:
            d14[idx - 1] += 1
        if event['n1'] == 0 and event['n2'] >= 1 and event['n3'] >= 1 \
            and event['n4'] == 0:
            d23[idx - 1] += 1
        if event['n1'] == 0 and event['n2'] >= 1 and event['n3'] == 0 \
            and event['n4'] >= 1:
            d24[idx - 1] += 1
        if event['n1'] == 0 and event['n2'] == 0 and event['n3'] >= 1 \
            and event['n4'] >= 1:
            d34[idx - 1] += 1

    figure()

    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    plot(x, pd2, label="Poisson")
    #plot(x, d12 / all, label="Data: 12")
    #plot(x, d13 / all, label="Data: 13")
    #plot(x, d14 / all, label="Data: 14")
    #plot(x, d23 / all, label="Data: 23")
    #plot(x, d24 / all, label="Data: 24")
    #plot(x, d34 / all, label="Data: 34")
    plot(x, (d12 + d13 + d14 + d23 + d24 + d34) / 6 / all,
         label="Simulation")

    xlabel("Mean density")
    ylabel("Probability")
    title("Prob: two detectors with particles")
    legend(loc='best')
    savefig("plots/poisson-two.pdf")

def plot_poisson_three():
    events = data.root.analysis.angle_0
  
    #k_events = kdata.getNode('/efficiency', 'events')
    #q = '(.9e15 <= k_energy) & (k_energy < 1.1e15) & (k_zenith < .35)'
    #kevents = k_events.readWhere(q)

    #global kd
    #k3 = []
    #kd = mean(kevents['k_dens_e'] + kevents['k_dens_mu'], 1)
    #for d0, d1 in progressbar(zip(bins[:-1], bins[1:])):
    #    k = kevents.compress('(d0 <= kd) & (kd < d1)')
    #    k = k['pulseheights'] 
    #UNFINISHED!

    bins = linspace(0, 10, 51)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']

    all = zeros(len(bins) - 1)
    d123 = zeros(len(bins) - 1)
    d124 = zeros(len(bins) - 1)
    d134 = zeros(len(bins) - 1)
    d234 = zeros(len(bins) - 1)
    for idx, event in progressbar(zip(bins.searchsorted(dens),
                                      events[:])):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        if event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] >= 1 \
            and event['n4'] == 0:
            d123[idx - 1] += 1
        if event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] == 0 \
            and event['n4'] >= 1:
            d124[idx - 1] += 1
        if event['n1'] >= 1 and event['n2'] == 0 and event['n3'] >= 1 \
            and event['n4'] >= 1:
            d134[idx - 1] += 1
        if event['n1'] == 0 and event['n2'] >= 1 and event['n3'] >= 1 \
            and event['n4'] >= 1:
            d234[idx - 1] += 1

    figure()

    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit

    plot(x, pd3, label="Poisson")
    #plot(x, d123 / all, label="Data: 123")
    #plot(x, d124 / all, label="Data: 124")
    #plot(x, d134 / all, label="Data: 134")
    #plot(x, d234 / all, label="Data: 234")
    plot(x, (d123 + d124 + d134 + d234) / 4 / all, label="Simulation")

    xlabel("Mean density")
    ylabel("Probability")
    title("Prob: three detectors with particles")
    legend(loc='best')
    savefig("plots/poisson-three.pdf")

def plot_poisson_four():
    events = data.root.analysis.angle_0

    global all, d1234, x
   
    bins = linspace(0, 10, 51)
    x = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])

    dens = events[:]['d1']

    all = zeros(len(bins) - 1)
    d1234 = zeros(len(bins) - 1)
    for idx, event in progressbar(zip(bins.searchsorted(dens),
                                      events[:])):
        try:
            all[idx - 1] += 1
        except IndexError:
            continue

        if event['n1'] >= 1 and event['n2'] >= 1 and event['n3'] >= 1 \
            and event['n4'] >= 1:
            d1234[idx - 1] += 1

    figure()

    global p0, pp, pd0, pd1, pd2, pd3, pd4

    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    #pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit
    pd4 = pp ** 4             # 4 detectors hit

    plot(x, pd4, label="Poisson")
    plot(x, d1234 / all, label="Simulation")

    xlabel("Mean density")
    ylabel("Probability")
    title("Prob: four detectors with particles")
    legend(loc='best')
    savefig("plots/poisson-four.pdf")


if __name__ == '__main__':
    try:
        data, kdata
    except NameError:
        data = tables.openFile(DATAFILE, 'r')
        kdata = tables.openFile(KASCADEDATAFILE, 'r')

    main()
