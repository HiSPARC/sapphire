import datetime
import urllib2
import json
import StringIO
import warnings

import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit


API_URL = 'http://data.hisparc.nl/api/stations/data/%d/%d/%d/'
HIST_URL = 'http://data.hisparc.nl/show/source/pulseintegral/%d/%d/%d/%d/'

MPV_FIT_WIDTH_FACTOR = .4


def main():
    today = datetime.date.today()
    yesterday = today - datetime.timedelta(days=1)
    station_ids = get_station_ids_with_data(yesterday)

    for station in station_ids[:1]:
        if station == 10:
            continue
        print station
        n, bins = get_histogram_for_station_on_date(station, yesterday)
        find_mpv = FindMostProbableValue()
        mpv, is_fitted = find_mpv.find_mpv_in_histogram(n, bins)

        figure()
        plot((bins[:-1] + bins[1:]) / 2., n)
        if is_fitted:
            axvline(mpv, c='g')
        else:
            axvline(mpv, c='r')
        title(station)


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


class FindMostProbableValue:
    def find_mpv_in_histogram(self, n, bins):
        first_guess = self.find_first_guess_mpv_in_histogram(n, bins)
        try:
            mpv = self.fit_mpv_in_histogram(n, bins, first_guess)
        except RuntimeError:
            warnings.warn("Fit failed, using first guess")
            return first_guess, False
        else:
            return mpv, True

    def find_first_guess_mpv_in_histogram(self, n, bins):
        """First guesst of most probable value in histogram.

        Algorithm: First: from the left: find the greatest value and
        cut off all data to the left of that maximum.  Now, you've cut off the
        where the trigger cuts in data.  Work from the right: find the
        location with the greatest decrease.  Then find the location of the
        maximum to the right of this location.

        """
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

    def fit_mpv_in_histogram(self, n, bins, first_guess):
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


if __name__ == '__main__':
    main()

#    today = datetime.date.today()
#    yesterday = today - datetime.timedelta(days=1)
#    figure()
#    n, bins = get_histogram_for_station_on_date(501, yesterday)
#    plot((bins[:-1] + bins[1:]) / 2., n)
#    mpv = find_first_guess_mpv_in_histogram(n, bins)
#    fit_mpv = fit_mpv_in_histogram(n, bins, mpv)
#    print mpv, fit_mpv
#    axvline(mpv, c='r')
#    axvline(fit_mpv, c='g')
