"""Find the most probable value in a HiSPARC spectrum.

:class:`FastFindMostProbableValueInSpectrum`
    find the most probable value in a HiSPARC spectrum using a fast
    algorithm
:class:`FindMostProbableValueInSpectrum`
    find the most probable value in a HiSPARC spectrum using a more
    elaborate approach

"""
import datetime
import urllib2
import json
import StringIO
import warnings

import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
import pylab as plt

from sapphire.api import Network


HIST_URL = 'http://data.hisparc.nl/show/source/pulseintegral/%d/%d/%d/%d/'

MPV_FIT_WIDTH_FACTOR = .4


class FastFindMostProbableValueInSpectrum(object):

    """Find the most probable value (MPV) in a HiSPARC spectrum.

    This is a fast algorithm to find the MPV value in a HiSPARC spectrum.
    The MPV value indicates the position of the minimum-ionizing particles
    (MIP) peak.  The algorithm makes some assumptions about the shape of
    the spectrum:

       * the spectrum includes the gamma peak (left-most part of
         spectrum) which has more counts per bin than the MIP peak.
       * ignoring the gamma peak, the MIP peak can be bracketed on the
         left by the numerically largest bin-to-bin increase in the
         number of counts.
       * the MIP peak can be approximated by a normal distribution.

    This algorithm can fail when the MIP peak is not pronounced.

    Public methods:

    :meth:`find_mpv`
        Find the most probable value
    :meth:`find_first_guess_mpv`
        Make a first guess of the most probable value
    :meth:`fit_mpv`
        Based on a first guess, fit the MIP peak to obtain the MPV

    """

    def __init__(self, n, bins):
        """Initialize the class instance.

        :param n, bins: histogram counts and bins, as obtained using
            :func:`numpy.histogram`.

        """
        self.n, self.bins = n, bins

    def find_mpv(self):
        """Find the most probable value.

        First perform a first guess, then use that value to fit the MIP
        peak.

        :return mpv: best guess of the most probable value
        :return boolean is_fitted: indicates if the fit was successful.

        """
        first_guess = self.find_first_guess_mpv()
        try:
            mpv = self.fit_mpv(first_guess)
        except RuntimeError:
            warnings.warn("Fit failed")
            return -999, False
        else:
            return mpv, True

    def find_first_guess_mpv(self):
        """First guess of most probable value.

        The algorithm is fast and simple. The following steps are
        performed:

           * From the left: find the greatest value and cut off all data
             to the left of that maximum.  We now assume the first
             datapoint to be the maximum of the gamma peak.
           * From the right: find the location of the greatest decrease
             from bin to bin.  We assume that this value is where the MIP
             peak dips before joining the gamma peak.
           * Find the maximum *to the right* of this value.  We assume
             this to be the approximate location of the MIP peak.

        :returns mpv: first guess of the most probable value

        """
        n, bins = self.n, self.bins

        # cut off trigger from the left
        left_idx = n.argmax()
        cut_n = n[left_idx:]

        # find greatest decrease from right
        delta_n = cut_n[:-1] - cut_n[1:]
        idx_greatest_decrease = delta_n.argmin()
        cut_cut_n = cut_n[idx_greatest_decrease:]

        # estimate most probable value with maximum in bracketed data
        idx_right_max = cut_cut_n.argmax()

        # calculate position of most probable value
        idx_mpv = idx_right_max + idx_greatest_decrease + left_idx
        mpv = (bins[idx_mpv] + bins[idx_mpv + 1]) / 2.

        return mpv

    def fit_mpv(self, first_guess, width_factor=MPV_FIT_WIDTH_FACTOR):
        """Fit a normal distribution to the MIP peak to obtain the MPV.

        A normal distribution is fitted to the spectrum in a restricted
        domain around the first guess value.  The width of the domain can
        be adjusted by the width_factor parameter.

        :param first_guess: approximate location of the most probable
            value
        :param width_factor: float in the range [0., 1.] to indicate the
            width of the fit domain.  The domain is given by
            [(1. - width_factor) * first_guess, (1. + width_factor) *
            first_guess]
        :returns mpv: mpv value obtained from the fit

        """
        n, bins = self.n, self.bins

        bins_x = (bins[:-1] + bins[1:]) / 2.

        # calculate fit domain
        left = (1. - width_factor) * first_guess
        right = (1. + width_factor) * first_guess

        # bracket histogram data
        x = bins_x.compress((left <= bins_x) & (bins_x < right))
        y = n.compress((left <= bins_x) & (bins_x < right))

        # sanity check: number of data points must be at least equal to
        # the number of fit parameters
        if len(x) < 3:
            raise RuntimeError("Number of data points not sufficient")

        # fit to a normal distribution
        f = lambda x, N, a, b: N * norm.pdf(x, loc=a, scale=b)
        popt, pcov = curve_fit(f, x, y, p0=(y.max(), first_guess,
                                            first_guess))
        mpv = popt[1]

        # sanity check: if MPV is outside domain, the MIP peak was not
        # bracketed correctly
        if mpv < x[0] or mpv > x[-1]:
            raise RuntimeError("Fitted MPV value outside range")

        return mpv


class FindMostProbableValueInSpectrum(FastFindMostProbableValueInSpectrum):

    """Find the most probable value (MPV) in a HiSPARC spectrum.

    This is an improved algorithm to find the MPV value in a HiSPARC spectrum.
    The MPV value indicates the position of the minimum-ionizing particles
    (MIP) peak.  The algorithm makes some assumptions about the shape of
    the spectrum:

       * the spectrum includes the gamma peak (left-most part of
         spectrum) which has more counts per bin than the MIP peak.
       * the spectrum can be fit with a power law, and the result loosely
         follows the spectrum over the complete range
       * when the power law fit is subtracted from the spectrum, the MIP
         peak either shows up clearly, or very poorly (when absent).

    This algorithm detects MIP peaks even if there is no clear dip between
    the gamma slope and MIP peak.  It even succeeds for most 'I think that
    might be a bump there' situations.

    Public methods:

    :meth:`find_mpv`
       Find the most probable value
    :meth:`find_first_guess_mpv`
       Make a first guess of the most probable value
    :meth:`fit_mpv`
       Based on a first guess, fit the MIP peak to obtain the MPV

    """

    def __init__(self, n, bins):
        """Initialize the class instance.

        On initialization, a power law is fitted to and subtracted from
        the data.  Then, the data is normalized so its maximum is 1.

        :param n, bins: histogram counts and bins, as obtained using
            :func:`numpy.histogram`.

        """
        reduced_y = self.fit_and_subtract_power_law(n, bins)
        super(FindMostProbableValueInSpectrum, self).__init__(reduced_y,
                                                              bins)

    def fit_and_subtract_power_law(self, n, bins):
        x = (bins[:-1] + bins[1:]) / 2

        # cut off trigger from the left
        left_idx = n.argmax()
        y = n[left_idx:]
        x = x[left_idx:]

        # Fit with a power law
        f = lambda x, a, b: a * x ** b
        popt, pcov = curve_fit(f, x, y)

        # Subtract power law from data
        ty = y - f(x, *popt)
        reduced_y = ty.clip(0, max(ty))

        # Fake single high bin at left, because the 'fast' algorithm will
        # cut off the maximum
        reduced_y[0] = max(reduced_y) + 1

        return reduced_y / reduced_y.max()

    def find_mpv(self):
        """Find the most probable value.

        First perform a first guess, then use that value to fit the MIP
        peak.

        :return mpv: best guess of the most probable value
        :return boolean is_fitted: indicates if the fit was successful.

        """
        smoothness = self.get_smoothness_of_data()
        if smoothness > 1:
            # Do not even try to fit, there is no MIP peak.
            return (-999, False)
        else:
            return super(FindMostProbableValueInSpectrum, self).find_mpv()

    def get_smoothness_of_data(self):
        """Quantify the smoothness of the data.

        After the power law fit and subtraction, the resulting data should
        either show a smooth bump (the MIP peak) or have a very ragged
        appearance.  This method quantifies the smoothness of the data.

        A smoothness < 1 is good.  Larger than 1 is bad.

        """
        # Do not include the fake maximum value in the first bin
        y = self.n[1:]
        # Based on two values, linearly extrapolate the third value.
        # Calculate the difference with the actual value and sum over
        # them.
        return ((2 * y[1:-1] - y[:-2] - y[2:]) ** 2).sum()


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
       :func:`numpy.histogram`.

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
