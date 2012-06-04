import datetime
import calendar

import tables
from pylab import *
from scipy.optimize import curve_fit

from sapphire.kascade import KascadeCoincidences
from artist import GraphArtist

import utils


USE_TEX = True

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


def do_matching_plots(data):
    #plot_nearest_neighbors(data)
    plot_residual_time_differences(data)

def plot_nearest_neighbors(data, limit=None):
    global coincidences
    hisparc_group = data.root.hisparc.cluster_kascade.station_601
    kascade_group = data.root.kascade

    coincidences = KascadeCoincidences(data, hisparc_group, kascade_group,
                                       ignore_existing=True)

    dt_opt = find_optimum_dt(coincidences, p0=-13, limit=1000)
    print dt_opt

    uncorrelated = None
    figure()
    for shift in -12, -13, dt_opt, -14:
        print "Shifting", shift
        coincidences.search_coincidences(shift, dtlimit=1, limit=limit)
        print "."
        dts = coincidences.coincidences['dt']
        n, bins, p = hist(abs(dts) / 1e9, bins=100, histtype='step',
                          label='%.3f s' % shift)
        if uncorrelated is None:
            uncorrelated = n, bins

    y, bins = uncorrelated
    x = (bins[:-1] + bins[1:]) / 2
    f = lambda x, N, a: N * exp(-a * x)
    popt, pcov = curve_fit(f, x, y)
    plot(x, f(x, *popt), label=r"$\lambda = %.2f$ Hz" % popt[1])

    yscale('log')
    xlabel("Time difference [s]")
    ylabel("Counts")
    legend()
    utils.saveplot()

def find_optimum_dt(coincidences, p0, delta=1e-3, limit=None):
    dt_new = p0
    dt_old = Inf

    while abs(dt_new - dt_old) > delta:
        print "Trying:", dt_new
        coincidences.search_coincidences(dt_new, dtlimit=1, limit=limit)
        dts = coincidences.coincidences['dt'] / 1e9

        dt_old = dt_new
        if median(dts) > 0:
            dt_new -= median(dts.compress(dts > 0))
        else:
            dt_new -= median(dts.compress(dts < 0))

    return dt_new

def plot_residual_time_differences(data):
    global idxes, dts
    events = data.root.kascade.events
    c_index = data.root.kascade.c_index

    t0 = make_timestamp(2008, 7, 2)
    t1 = make_timestamp(2008, 7, 3)

    idxes = events.getWhereList('(t0 <= timestamp) & (timestamp < t1)')
    t0_idx = min(idxes)
    t1_idx = max(idxes)

    dts = c_index.readWhere('(t0_idx <= k_idx) & (k_idx < t1_idx)',
                            field='dt')
    all_dts = c_index.col('dt')

    figure()
    subplot(121)
    hist(all_dts / 1e3, bins=arange(-10, 2, .01), histtype='step')
    title("July 1 - Aug 6, 2008")
    xlabel("Time difference [us]")
    ylabel("Counts")

    subplot(122)
    hist(dts / 1e3, bins=arange(-8, -6, .01), histtype='step')
    title("July 2, 2008")
    xlabel("Time difference [us]")
    utils.saveplot()

    graph = GraphArtist(width=r'.45\linewidth')
    n, bins = histogram(all_dts / 1e3, bins=arange(-10, 2, .01))
    graph.histogram(n, bins)
    graph.set_title("Jul 1 - Aug 6, 2008")
    graph.set_xlabel(r"Time difference [\si{\micro\second}]")
    graph.set_ylabel("Counts")
    graph.set_xlimits(-10, 2)
    graph.set_ylimits(min=0)
    graph.save('plots/MAT-residual-time-differences-weeks')

    graph = GraphArtist(width=r'.45\linewidth')
    n, bins = histogram(dts / 1e3, bins=arange(-8, -6, .01))
    graph.histogram(n, bins)
    graph.set_title("Jul 2, 2008")
    graph.set_xlabel(r"Time difference [\si{\micro\second}]")
    graph.set_ylabel("Counts")
    graph.set_xlimits(-8, -6)
    graph.set_ylimits(min=0)
    graph.save('plots/MAT-residual-time-differences-day')

def make_timestamp(year, month, day, hour=0, minutes=0, seconds=0):
    dt = datetime.datetime(year, month, day, hour, minutes, seconds)
    timetuple = dt.utctimetuple()
    timestamp = calendar.timegm(timetuple)
    return timestamp


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'r')

    utils.set_prefix('MAT-')
    do_matching_plots(data)
