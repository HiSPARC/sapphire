import datetime
import urllib2
import json
import StringIO
import warnings

import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit

import pylab as plt


API_URL = 'http://data.hisparc.nl/api/stations/data/%d/%d/%d/'
HIST_URL = 'http://data.hisparc.nl/show/source/pulseintegral/%d/%d/%d/%d/'

MPV_FIT_WIDTH_FACTOR = .4


class FindMostProbableValue:
    def __init__(self, n, bins):
        self.n, self.bins = n, bins

    def find_mpv_in_histogram(self):
        first_guess = self.find_first_guess_mpv_in_histogram()
        try:
            mpv = self.fit_mpv_in_histogram(first_guess)
        except RuntimeError:
            warnings.warn("Fit failed, using first guess")
            return first_guess, False
        else:
            return mpv, True

    def find_first_guess_mpv_in_histogram(self):
        """First guesst of most probable value in histogram.

        Algorithm: First: from the left: find the greatest value and
        cut off all data to the left of that maximum.  Now, you've cut off the
        where the trigger cuts in data.  Work from the right: find the
        location with the greatest decrease.  Then find the location of the
        maximum to the right of this location.

        """
        n, bins = self.n, self.bins

        # cut off trigger from the left
        left_idx = n.argmax()
        cut_n = n[left_idx:]

        # find mpv peak from right
        delta_n = cut_n[:-1] - cut_n[1:]
        idx_greatest_decrease = delta_n.argmin()
        idx_right_max = cut_n[idx_greatest_decrease:].argmax()

        idx_mpv = idx_right_max + idx_greatest_decrease + left_idx
        mpv = (bins[idx_mpv] + bins[idx_mpv + 1]) / 2.

        return mpv

    def fit_mpv_in_histogram(self, first_guess):
        n, bins = self.n, self.bins

        bins_x = (bins[:-1] + bins[1:]) / 2.

        left = (1. - MPV_FIT_WIDTH_FACTOR) * first_guess
        right = (1. + MPV_FIT_WIDTH_FACTOR) * first_guess

        x = bins_x.compress((left <= bins_x) & (bins_x < right))
        y = n.compress((left <= bins_x) & (bins_x < right))

        f = lambda x, N, a, b: N * norm.pdf(x, loc=a, scale=b)
        popt, pcov = curve_fit(f, x, y, p0=(y.max(), first_guess,
                                            first_guess))
        mpv = popt[1]

        if mpv < x[0] or mpv > x[-1]:
            raise RuntimeError("Fitted MPV value outside range")

        return mpv


def main():
    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    station_ids = get_station_ids_with_data(yesterday)

    for station in station_ids:
        if station == 10:
            continue
        print station
        n, bins = get_histogram_for_station_on_date(station, yesterday)
        find_mpv = FindMostProbableValue(n, bins)
        mpv, is_fitted = find_mpv.find_mpv_in_histogram()

        plt.figure()
        plt.plot((bins[:-1] + bins[1:]) / 2., n)
        if is_fitted:
            plt.axvline(mpv, c='g')
        else:
            plt.axvline(mpv, c='r')
        plt.title(station)


def get_station_ids_with_data(date):
    url = API_URL % (date.year, date.month, date.day)

    reply = urllib2.urlopen(url)
    reply = reply.read()

    station_list = json.loads(reply)
    station_ids = [int(u['number']) for u in station_list]

    return station_ids


def get_histogram_for_station_on_date(station_id, date):
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
