import tables
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import scoreatpercentile
from pylab import *

import utils

from artist import GraphArtist


USE_TEX = False

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


def main():
    global data
    data = tables.open_file('master-ch4v2.h5', 'r')
    #utils.set_suffix('E_1PeV')

    #scatterplot_core_distance_vs_time()
    #median_core_distance_vs_time()
    boxplot_core_distance_vs_time()
    #hists_core_distance_vs_time()
    plot_front_passage()


def scatterplot_core_distance_vs_time():
    plt.figure()

    sim = data.root.showers.E_1PeV.zenith_0
    electrons = sim.electrons

    plt.loglog(electrons[:]['core_distance'], electrons[:]['arrival_time'], ',')
    plt.xlim(1e0, 1e2)
    plt.ylim(1e-3, 1e3)

    plt.xlabel("Core distance [m]")
    plt.ylabel("Arrival time [ns]")

    utils.title("Shower front timing structure")
    utils.saveplot()


def median_core_distance_vs_time():
    plt.figure()
    plot_and_fit_statistic(lambda a: scoreatpercentile(a, 25))
    plot_and_fit_statistic(lambda a: scoreatpercentile(a, 75))

    utils.title("Shower front timing structure (25, 75 %)")
    utils.saveplot()
    plt.xlabel("Core distance [m]")
    plt.ylabel("Median arrival time [ns]")
    legend(loc='lower right')


def plot_and_fit_statistic(func):
    sim = data.root.showers.E_1PeV.zenith_0
    electrons = sim.electrons

    bins = np.logspace(0, 2, 25)
    x, y = [], []
    for low, high in zip(bins[:-1], bins[1:]):
        sel = electrons.read_where('(low < core_distance) & (core_distance <= high)')
        statistic = func(sel[:]['arrival_time'])
        x.append(np.mean([low, high]))
        y.append(statistic)

    plt.loglog(x, y)

    logx = log10(x)
    logy = log10(y)
    logf = lambda x, a, b: a * x + b
    g = lambda x, a, b: 10 ** logf(log10(x), a, b)
    popt, pcov = curve_fit(logf, logx, logy)
    plot(x, g(x, *popt), label="f(x) = %.2e * x ^ %.2e" % (10 ** popt[1],
                                                           popt[0]))


def boxplot_core_distance_vs_time():
    plt.figure()

    sim = data.root.showers.E_1PeV.zenith_0.shower_0
    leptons = sim.leptons

    #bins = np.logspace(0, 2, 25)
    bins = np.linspace(0, 100, 15)
    x, arrival_time, widths = [], [], []
    t25, t50, t75 = [], [], []
    for low, high in zip(bins[:-1], bins[1:]):
        sel = leptons.read_where('(low < core_distance) & (core_distance <= high)')
        x.append(np.mean([low, high]))
        arrival_time.append(sel[:]['arrival_time'])
        widths.append((high - low) / 2)
        ts = sel[:]['arrival_time']
        t25.append(scoreatpercentile(ts, 25))
        t50.append(scoreatpercentile(ts, 50))
        t75.append(scoreatpercentile(ts, 75))

    fill_between(x, t25, t75, color='0.75')
    plot(x, t50, 'o-', color='black')

    plt.xlabel("Core distance [m]")
    plt.ylabel("Arrival time [ns]")

    #utils.title("Shower front timing structure")
    utils.saveplot()

    graph = GraphArtist()
    graph.plot(x, t50, linestyle=None)
    graph.shade_region(x, t25, t75)
    graph.set_xlabel(r"Core distance [\si{\meter}]")
    graph.set_ylabel(r"Arrival time [\si{\nano\second}]")
    graph.set_ylimits(0, 30)
    graph.set_xlimits(0, 100)
    graph.save('plots/front-passage-vs-R')


def hists_core_distance_vs_time():
    plt.figure()

    sim = data.root.showers.E_1PeV.zenith_0
    electrons = sim.electrons

    bins = np.logspace(0, 2, 5)
    for low, high in zip(bins[:-1], bins[1:]):
        sel = electrons.read_where('(low < core_distance) & (core_distance <= high)')
        arrival_time = sel[:]['arrival_time']
        plt.hist(arrival_time, bins=np.logspace(-2, 3, 50), histtype='step',
                 label="%.2f <= log10(R) < %.2f" % (np.log10(low),
                                                    np.log10(high)))

    plt.xscale('log')

    plt.xlabel("Arrival Time [ns]")
    plt.ylabel("Count")
    plt.legend(loc='upper left')

    utils.title("Shower front timing structure")
    utils.saveplot()


def plot_front_passage():
    sim = data.root.showers.E_1PeV.zenith_0.shower_0
    leptons = sim.leptons
    R = 40
    dR = 2
    low = R - dR
    high = R + dR
    global t
    t = leptons.read_where('(low < core_distance) & (core_distance <= high)',
                          field='arrival_time')

    n, bins, patches = hist(t, bins=linspace(0, 30, 31), histtype='step')

    graph = GraphArtist()
    graph.histogram(n, bins)
    graph.set_xlabel(r"Arrival time [\si{\nano\second}]")
    graph.set_ylabel("Number of leptons")
    graph.set_ylimits(min=0)
    graph.set_xlimits(0, 30)
    graph.save('plots/front-passage')


if __name__ == '__main__':
    main()
