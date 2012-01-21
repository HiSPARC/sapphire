from __future__ import division

import tables
import numpy as np
import pylab as plt
from scipy import optimize

from sapphire.analysis import landau


RANGE_MAX = 40000
N_BINS = 400


class ReconstructionEfficiency(object):
    def __init__(self, data):
        self.data = data
        self.scintillator = landau.Scintillator()

    def main(self):
        self.plot_landau_fit()

    def calc_charged_fraction(self, x, y, p_gamma, p_landau):
        y_reduced = y - self.gamma_func(x, *p_gamma)

        mev_scale = p_landau[1]
        max_pos = 3.38 / mev_scale
        print "max_pos", max_pos

        y_landau = self.scintillator.conv_landau_for_x(x, *p_landau)
        y_charged_left = y_landau.compress(x <= max_pos)
        y_charged_right = y_reduced.compress(max_pos < x)
        y_charged = array(y_charged_left.tolist() +
                          y_charged_right.tolist())

        plt.plot(x, y_charged, '-')

        N_full = y.sum()
        N_charged = y_charged.sum()

        print "full, charged, fraction:", N_full, N_charged, N_full / N_charged
        return N_charged / N_full

    def full_spectrum_fit(self, x, y, p0_gamma, p0_landau):
        p_gamma = self.fit_gammas_to_data(x, y, p0_gamma)
        p_landau = self.fit_conv_landau_to_data(x, y - self.gamma_func(x, *p_gamma),
                                                p0_landau)
        p_gamma, p_landau = self.fit_complete(x, y, p_gamma, p_landau)
        return p_gamma, p_landau

    def constrained_full_spectrum_fit(self, x, y, p0_gamma, p0_landau):
        p_gamma, p_landau = self.constrained_fit_complete(x, y, p0_gamma, p0_landau)
        return p_gamma, p_landau

    def plot_landau_fit(self):
        global x, n, bins, p_gamma, p_landau

        events = self.data.root.hisparc.cluster_kascade.station_601.events
        ph0 = events.col('integrals')[:, 0]

        bins = np.linspace(0, RANGE_MAX, N_BINS + 1)
        n, bins = np.histogram(ph0, bins=bins)
        x = (bins[:-1] + bins[1:]) / 2
    
        p_gamma, p_landau = self.full_spectrum_fit(x, n, (1., 1.),
                                                   (5e3 / .32, 3.38 / 5000, 1.))
        print "FULL FIT"
        print p_gamma, p_landau

        n /= 10
        p_gamma, p_landau = self.constrained_full_spectrum_fit(x, n, p_gamma, p_landau)
        print "CONSTRAINED FIT"
        print p_gamma, p_landau

        clf()
        print self.calc_charged_fraction(x, n, p_gamma, p_landau)

        plt.plot(x, n)
        self.plot_landau_and_gamma(x, p_gamma, p_landau)
        #plt.plot(x, n - self.gamma_func(x, *p_gamma))
        plt.yscale('log')
        plt.xlim(xmin=0)
        plt.ylim(ymin=1e1)

    def plot_landau_and_gamma(self, x, p_gamma, p_landau):
        gammas = self.gamma_func(x, *p_gamma)
        plt.plot(x, gammas)

        nx = np.linspace(-RANGE_MAX, RANGE_MAX, N_BINS * 2 + 1)
        nlandaus = self.scintillator.conv_landau(nx, *p_landau)
        landaus = np.interp(x, nx, nlandaus)
        plt.plot(x, landaus)

        plt.plot(x, gammas + landaus)


    def fit_gammas_to_data(self, x, y, p0):
        condition = (500 <= x) & (x < 2000)
        x_trunc = x.compress(condition)
        y_trunc = y.compress(condition)
        popt, pcov = optimize.curve_fit(self.gamma_func, x_trunc, y_trunc, p0=p0)
        return popt

    def gamma_func(self, x, N, a):
        return N * x ** -a

    def fit_conv_landau_to_data(self, x, y, p0):
        x_symm = np.linspace(-RANGE_MAX, RANGE_MAX, N_BINS * 2 + 1)
        y_symm = np.interp(x_symm, x, y)
        popt = optimize.fmin(self.scintillator.residuals, p0,
                             (x_symm, y_symm, 4500, 5500))
        return popt

    def fit_complete(self, x, y, p_gamma, p_landau):
        x_symm = np.linspace(-RANGE_MAX, RANGE_MAX, N_BINS * 2 + 1)
        y_symm = np.interp(x_symm, x, y)
        p0 = list(p_gamma) + list(p_landau)
        popt = optimize.fmin(self.complete_residuals, p0,
                             (self.scintillator, x_symm, y_symm, 500, 6000),
                             maxfun=100000)
        return popt[:2], popt[2:]

    def constrained_fit_complete(self, x, y, p_gamma, p_landau):
        x_symm = np.linspace(-RANGE_MAX, RANGE_MAX, N_BINS * 2 + 1)
        y_symm = np.interp(x_symm, x, y)
        N_gamma = p_gamma[0]
        N_landau = p_landau[0]
        popt = optimize.fmin(self.constrained_complete_residuals,
                             (N_gamma, N_landau),
                             (self.scintillator, x_symm, y_symm, p_gamma,
                              p_landau, 500, 6000),
                             maxfun=100000)
        p_gamma[0] = popt[0]
        p_landau[0] = popt[1]
        return p_gamma, p_landau

    def complete_residuals(self, par, scintillator, x, y, a, b):
        landaus = scintillator.conv_landau(x, *par[2:])
        gammas = self.gamma_func(x, *par[:2])
        residuals = (y - (gammas + landaus)) ** 2
        residuals = residuals.compress((a <= x) & (x < b))
        residuals = residuals.sum()
        return residuals

    def constrained_complete_residuals(self, par, scintillator, x, y,
                                       p_gamma, p_landau, a, b):
        full_par = (par[0], p_gamma[1], par[1], p_landau[1], p_landau[2])
        return self.complete_residuals(full_par, scintillator, x, y, a, b)


if __name__ == '__main__':
    np.seterr(invalid='ignore', divide='ignore')

    if 'data' not in globals():
        data = tables.openFile('kascade.h5', 'r')

    efficiency = ReconstructionEfficiency(data)
    efficiency.main()
