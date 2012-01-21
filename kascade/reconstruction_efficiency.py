import tables
import numpy as np
import pylab as plt
from scipy import optimize

from sapphire.analysis import landau


RANGE_MAX = 40000
N_BINS = 400


def main():
    plot_landau_fit()

def plot_landau_fit():
    global x, n, bins, p_gamma, p_landau

    events = data.root.hisparc.cluster_kascade.station_601.events
    ph0 = events.col('integrals')[:, 0]

    bins = np.linspace(0, RANGE_MAX, N_BINS + 1)
    n, bins = np.histogram(ph0, bins=bins)
    x = (bins[:-1] + bins[1:]) / 2

    p_gamma = fit_gammas_to_data(x, n)
    p_landau = fit_conv_landau_to_data(x, n - gamma_func(x, *p_gamma))
    p_gamma, p_landau = fit_complete(x, n, p_gamma, p_landau)

    clf()
    plt.plot(x, n)
    plot_landau_and_gamma(x, p_gamma, p_landau)
    plt.plot(x, n - gamma_func(x, *p_gamma))
    plt.yscale('log')
    plt.xlim(xmin=0)
    plt.ylim(ymin=1e2)

def plot_landau_and_gamma(x, p_gamma, p_landau):
    gammas = gamma_func(x, *p_gamma)
    plot(x, gammas)

    nx = linspace(-RANGE_MAX, RANGE_MAX, N_BINS * 2 + 1)
    nlandaus = scintillator.conv_landau(nx, *p_landau)
    landaus = interp(x, nx, nlandaus)
    plot(x, landaus)

    plot(x, gammas + landaus)

def fit_gammas_to_data(x, y):
    condition = (500 <= x) & (x < 2000)
    x_trunc = x.compress(condition)
    y_trunc = y.compress(condition)
    popt, pcov = optimize.curve_fit(gamma_func, x_trunc, y_trunc)
    return popt

gamma_func = lambda x, N, a: N * x ** -a

def fit_conv_landau_to_data(x, y, p0=(1e4 / .32, 3.38 / 5000, 1)):
    x_symm = np.linspace(-RANGE_MAX, RANGE_MAX, N_BINS * 2 + 1)
    y_symm = np.interp(x_symm, x, y)
    popt = optimize.fmin(scintillator.residuals, p0,
                         (x_symm, y_symm, 4500, 5500))
    return popt

def fit_complete(x, y, p_gamma, p_landau):
    x_symm = np.linspace(-RANGE_MAX, RANGE_MAX, N_BINS * 2 + 1)
    y_symm = np.interp(x_symm, x, y)
    p0 = list(p_gamma) + list(p_landau)
    popt = optimize.fmin(complete_residuals, p0,
                         (scintillator, x_symm, y_symm, 500, 6000),
                         maxfun=100000)
    print p0
    print popt, '%e' % complete_residuals(popt, scintillator, x_symm, y_symm, 500, 6000)
    print [u / v for u, v in zip(p0, popt)]
    bounds = [(0.7 * u, 1.3 * u) for u in p0]
    popt, n, r = optimize.fmin_tnc(complete_residuals, p0, bounds=bounds,
                                        approx_grad=True, maxfun=10000,
                                        args=(scintillator,
                                              x_symm, y_symm, 500, 6000))
    print popt, '%e' % complete_residuals(popt, scintillator, x_symm, y_symm, 500, 6000)
    print n

    print
    print
    print [u / v for u, v in zip(p0, popt)]
    return popt[:2], popt[2:]

def complete_residuals(par, scintillator, x, y, a, b):
    landaus = scintillator.conv_landau(x, *par[2:])
    gammas = gamma_func(x, *par[:2])
    residuals = (y - (gammas + landaus)) ** 2
    residuals = residuals.compress((a <= x) & (x < b))
    residuals = residuals.sum()
    return residuals


if __name__ == '__main__':
    np.seterr(invalid='ignore', divide='ignore')

    if 'data' not in globals():
        data = tables.openFile('kascade.h5', 'r')
    if 'scintillator' not in globals():
        scintillator = landau.Scintillator()

    main()
