import tables
from scipy.optimize import curve_fit

from numpy import *
from pylab import *


def frac_bins(low, high, binsize, nbins=1):
    binsize = binsize * nbins
    low = low - .5 * binsize
    high = ceil((high - low) / binsize) * binsize + low + .5 * binsize
    return arange(low, high, binsize)


def fit_gauss_to_timings(events, timing_data):
    events = events[:]
    gauss = lambda x, N, mu, sigma: N * normpdf(x, mu, sigma)

    figure()
    for i, j in [(1, 0), (1, 2), (1, 3)]:
        dt = []
        for t, c in zip(timing_data, events):
            ph = c['pulseheights']
            if min([ph[i], ph[j]]) / 350. >= 2.:
                dt.append(t[i] - t[j])
        n, bins, patches = hist(dt, bins=frac_bins(-40, 40, 2.5),
                                histtype='step', label="detector %d - %d"
                                % (i + 1, j + 1))
        b = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
        popt, pcov = curve_fit(gauss, b, n)
        x = linspace(-40, 40, 100)
        plot(x, gauss(x, *popt), label="mu: %.2f, sigma: %.2f" % (popt[1],
                                                                  popt[2]))
    legend(prop={'size': 'small'})
    title("Delta t's for D >= 2.")
    xlabel("Time (ns)")
    ylabel("Count")


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.open_file('kascade.h5', 'r')

    events = data.root.kascade_new.coincidences
    timings = data.root.analysis_new.timing_data.read()
    timings_linear = data.root.analysis_new.timing_data_linear.read()

    fit_gauss_to_timings(events, timings)
    fit_gauss_to_timings(events, timings_linear)
