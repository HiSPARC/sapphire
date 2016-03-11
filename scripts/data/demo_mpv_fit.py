"""Demo the find_mpv module

Get the pulseintegral data from all stations with data yesterday.
Fit the MPV and show the results

"""
import datetime
import warnings

import numpy as np
import pylab as plt

from sapphire import Network, Station, FindMostProbableValueInSpectrum


COLORS = ['black', 'red', 'green', 'blue']


def main():
    """Demo the MPV finder with actual data."""

    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    station_ids = get_station_ids_with_data(yesterday)

    for station in station_ids:
        plt.figure()
        for did in range(Station(station).n_detectors()):
            n, bins = get_histogram_for_station_on_date(station, yesterday,
                                                        did)
            find_mpv = FindMostProbableValueInSpectrum(n, bins)
            mpv, is_fitted = find_mpv.find_mpv()

            plt.plot((bins[:-1] + bins[1:]) / 2., n, c=COLORS[did])
            lines = ['dotted', 'solid']
            plt.axvline(mpv + did * (bins[1] - bins[0]) / 20.,
                        c=COLORS[did], ls=lines[is_fitted])
            plt.title(station)
            plt.xlim(0, bins[len(bins)/2])
            plt.yscale('log')


def get_station_ids_with_data(date):
    """Return the station ids having data on this day."""

    station_list = Network.stations_with_data(date.year, date.month, date.day)
    station_ids = [int(station['number']) for station in station_list]

    return station_ids


def get_histogram_for_station_on_date(station_id, date, did):
    """Return a histogram of the spectrum of a station on a date.

    :return n, bins: histogram counts and bins, as obtained using
       ``numpy.histogram``.

    """
    data = Station(station_id).pulse_integral(date.year, date.month, date.day)

    bins = list(data['pulseintegral'])
    bins.append(bins[-1] + (bins[-1] - bins[-2]))
    bins = np.array(bins)

    n = data['pi%d' % (did + 1)]

    return n, bins


if __name__ == '__main__':
    warnings.simplefilter('always')
    main()
    plt.show()
