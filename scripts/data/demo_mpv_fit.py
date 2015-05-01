"""Demo the find_mpv module

Get the pulseintegral data from all stations with data yesterday.
Fit the MPV and show the results

"""
import datetime
import urllib2
import StringIO
import warnings

import numpy as np
import pylab as plt

from sapphire.analysis.find_mpv import FindMostProbableValueInSpectrum
from sapphire.api import Network


HIST_URL = 'http://data.hisparc.nl/show/source/pulseintegral/%d/%d/%d/%d/'


def main():
    """Demo the MPV finder with actual data."""

    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    station_ids = get_station_ids_with_data(yesterday)

    for station in station_ids:
        if station == 10:
            continue
        print station
        n, bins = get_histogram_for_station_on_date(station, yesterday)
        find_mpv = FindMostProbableValueInSpectrum(n, bins)
        mpv, is_fitted = find_mpv.find_mpv()

        plt.figure()
        plt.plot((bins[:-1] + bins[1:]) / 2., n)
        if is_fitted:
            plt.axvline(mpv, c='g')
        else:
            plt.axvline(mpv, c='r')
        plt.title(station)
        plt.yscale('log')


def get_station_ids_with_data(date):
    """Return the station ids having data on this day."""

    station_list = Network.stations_with_data(date.year, date.month, date.day)
    station_ids = [int(station['number']) for station in station_list]

    return station_ids


def get_histogram_for_station_on_date(station_id, date):
    """Return a histogram of the spectrum of a station on a date.

    :return n, bins: histogram counts and bins, as obtained using
       ``numpy.histogram``.

    """
    url = HIST_URL % (station_id, date.year, date.month, date.day)

    reply = urllib2.urlopen(url)
    reply = reply.read()

    file_like = StringIO.StringIO(reply)
    data = np.genfromtxt(file_like)

    bins = data[:, 0]
    bins = list(bins)
    bins.append(bins[-1] + (bins[-1] - bins[-2]))
    bins = np.array(bins)

    n = data[:, 1]

    return n, bins


if __name__ == '__main__':
    warnings.simplefilter('always')
    main()
    plt.show()
