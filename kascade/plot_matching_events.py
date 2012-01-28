import datetime
import calendar

import tables
from pylab import *

from sapphire.kascade import KascadeCoincidences

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
    plot_nearest_neighbors(data, limit=10000)
    plot_residual_time_differences(data)

def plot_nearest_neighbors(data, limit=None):
    global coincidences
    hisparc_group = data.root.hisparc.cluster_kascade.station_601
    kascade_group = data.root.kascade

    coincidences = KascadeCoincidences(data, hisparc_group, kascade_group,
                                       ignore_existing=True)

    dt_opt = find_optimum_dt(coincidences, p0=-13, limit=1000)
    print dt_opt

    figure()
    for shift in -12, -13, dt_opt, -14:
        coincidences.search_coincidences(shift, dtlimit=1, limit=limit)
        dts = coincidences.coincidences['dt']
        hist(abs(dts) / 1e9, bins=100, histtype='step', label='%.3f' % shift)

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

    t0 = make_timestamp(2008, 6, 4)
    t1 = make_timestamp(2008, 8, 4, 1)

    idxes = events.getWhereList('(t0 <= timestamp) & (timestamp < t1)')
    t0_idx = min(idxes)
    t1_idx = max(idxes)

    dts = c_index.readWhere('(t0_idx <= k_idx) & (k_idx < t1_idx)',
                            field='dt')

    figure()
    hist(dts / 1e3, bins=200, histtype='step')
    xlabel("Time difference [ns]")
    ylabel("Counts")
    utils.saveplot()

def make_timestamp(year, month, day, hour=0, minutes=0, seconds=0):
    dt = datetime.datetime(year, month, day, hour, minutes, seconds)
    timetuple = dt.utctimetuple()
    timestamp = calendar.timegm(timetuple)
    return timestamp


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade-small.h5', 'r')

    utils.set_prefix('MAT')
    do_matching_plots(data)
