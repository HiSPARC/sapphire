import tables
import scipy

from sapphire.analysis import landau


def main(data):
    s = landau.Scintillator()
    for num_detector in range(1, 5):
        fit_to_pulseheights(data, s, num_detector)
        fit_to_integrals(data, s, num_detector)
        fit_using_just_gauss_to_integrals(data, num_detector)
        print


def fit_to_pulseheights(data, s, num_detector):
    figure()

    events = data.root.hisparc.cluster_kascade.station_601.events

    ph = events.col('pulseheights')[:, num_detector - 1]
    print "Fitted to pulseheights, detector", num_detector
    popt = do_fit_to_data(ph, s, 1000, 201, (380, .3 * 380))
    print "Relative Gauss width (1st try):", popt[2] / 3.38
    center = 3.38 / popt[1]
    popt = do_fit_to_data(ph, s, 1000, 201, (center, .3 * center))
    print "Relative Gauss width (2nd try):", popt[2] / 3.38


def fit_to_integrals(data, s, num_detector):
    figure()

    events = data.root.hisparc.cluster_kascade.station_601.events

    intg = events.col('integrals')[:, num_detector - 1]
    print "Fitted to integrals, detector", num_detector
    popt = do_fit_to_data(intg, s, 20000, 201, (5000, .3 * 5000))
    print "Relative Gauss width (1st try):", popt[2] / 3.38
    center = 3.38 / popt[1]
    popt = do_fit_to_data(intg, s, 20000, 201, (center, .3 * center))
    print "Relative Gauss width (2nd try):", popt[2] / 3.38


def fit_using_just_gauss_to_integrals(data, num_detector):
    figure()

    events = data.root.hisparc.cluster_kascade.station_601.events

    intg = events.col('integrals')[:, num_detector - 1]
    print "Fitted to integrals (Gauss only), detector", num_detector
    popt = do_fit_to_data_using_gauss(intg, 20000, 201,
                                      (5000, .3 * 5000))
    print "Relative Gauss width (1st try):", popt[2] / 3.38
    center = 3.38 / popt[1]
    popt = do_fit_to_data_using_gauss(intg, 20000, 201,
                                      (center, .3 * center))
    print "Relative Gauss width (2nd try):", popt[2] / 3.38


def do_fit_to_data(dataset, s, max_hist_value, n_bins, guess):
    center, width = guess

    n, bins, patches = hist(dataset, bins=linspace(0, max_hist_value,
                                                   n_bins, 'b'),
                            histtype='step')
    yscale('log')

    x = (bins[:-1] + bins[1:]) / 2
    y = n
    low = center - width
    high = center + width
    sx = x.compress((low <= x) & (x < high))
    sy = y.compress((low <= x) & (x < high))

    guess_count = interp(center, sx, sy)

    popt, pcov = scipy.optimize.curve_fit(s.conv_landau_for_x, sx, sy,
                                          p0=(guess_count, 3.38 / center,
                                              1.))

    plot(sx, sy, 'r')
    plot(x, s.conv_landau_for_x(x, *popt), 'g')

    return popt


def do_fit_to_data_using_gauss(dataset, max_hist_value, n_bins, guess):
    center, width = guess

    n, bins, patches = hist(dataset, bins=linspace(0, max_hist_value,
                                                   n_bins, 'b'),
                            histtype='step')
    yscale('log')

    x = (bins[:-1] + bins[1:]) / 2
    y = n
    low = center - width
    high = center + width
    sx = x.compress((low <= x) & (x < high))
    sy = y.compress((low <= x) & (x < high))

    guess_count = interp(center, sx, sy)

    f = lambda u, N, scale, sigma: N * scipy.stats.norm.pdf(u * scale,
                                                         loc=3.38,
                                                         scale=sigma)
    popt, pcov = scipy.optimize.curve_fit(f, sx, sy,
                                          p0=(guess_count, 3.38 / center,
                                              1.))

    plot(sx, sy, 'r')
    plot(x, f(x, *popt), 'g')
    ylim(ymin=1e1)

    print popt

    return popt


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.open_file('kascade.h5')

    main(data)
