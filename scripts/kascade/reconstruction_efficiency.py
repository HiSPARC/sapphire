from __future__ import division

import tables
import numpy as np
import pylab as plt
from scipy import optimize, stats

from sapphire.analysis import landau

import utils

from artist import GraphArtist
import artist.utils


RANGE_MAX = 40000
N_BINS = 400

LOW, HIGH = 500, 5500

VNS = .57e-3 * 2.5

USE_TEX = True

# For matplotlib plots
if USE_TEX:
    plt.rcParams['font.serif'] = 'Computer Modern'
    plt.rcParams['font.sans-serif'] = 'Computer Modern'
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['figure.figsize'] = [4 * x for x in (1, 2. / 3)]
    plt.rcParams['figure.subplot.left'] = 0.175
    plt.rcParams['figure.subplot.bottom'] = 0.175
    plt.rcParams['font.size'] = 10
    plt.rcParams['legend.fontsize'] = 'small'
    plt.rcParams['text.usetex'] = True


class ReconstructionEfficiency(object):
    def __init__(self, data):
        global scintillator
        self.data = data

        if 'scintillator' in globals():
            self.scintillator = scintillator
        else:
            self.scintillator = landau.Scintillator()
            scintillator = self.scintillator

    def main(self):
        self.plot_spectrum_fit_chisq()
        self.plot_gamma_landau_fit()
        self.plot_detection_efficiency()

    def calc_charged_fraction(self, x, y, p_gamma, p_landau):
        y_charged = self.calc_charged_spectrum(x, y, p_gamma, p_landau)

        N_full = y.sum()
        N_charged = y_charged.sum()

        return N_charged / N_full

    def calc_charged_spectrum(self, x, y, p_gamma, p_landau):
        y_landau = self.scintillator.conv_landau_for_x(x, *p_landau)
        max_pos = x[y_landau.argmax()]

        y_gamma = self.gamma_func(x, *p_gamma)
        y_gamma_trunc = np.where(x <= 3 * max_pos, y_gamma, 0.)

        y_reduced = y - y_gamma_trunc

        mev_scale = p_landau[1]

        y_charged_left = y_landau.compress(x <= max_pos)
        y_charged_right = y_reduced.compress(max_pos < x)
        y_charged = np.array(y_charged_left.tolist() +
                             y_charged_right.tolist())

        return y_charged

    def full_spectrum_fit(self, x, y, p0_gamma, p0_landau):
        p_gamma = self.fit_gammas_to_data(x, y, p0_gamma)
        p_landau = self.fit_conv_landau_to_data(x, y - self.gamma_func(x, *p_gamma),
                                                p0_landau)
        p_gamma, p_landau = self.fit_complete(x, y, p_gamma, p_landau)
        return p_gamma, p_landau

    def constrained_full_spectrum_fit(self, x, y, p0_gamma, p0_landau):
        p_gamma, p_landau = self.constrained_fit_complete(x, y, p0_gamma, p0_landau)
        return p_gamma, p_landau

    def plot_gamma_landau_fit(self):
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

        plt.figure()
        print self.calc_charged_fraction(x, n, p_gamma, p_landau)

        plt.plot(x * VNS, n)
        self.plot_landau_and_gamma(x, p_gamma, p_landau)
        #plt.plot(x, n - self.gamma_func(x, *p_gamma))
        plt.xlabel("Pulse integral [V ns]")
        plt.ylabel("Count")
        plt.yscale('log')
        plt.xlim(0, 30)
        plt.ylim(1e1, 1e4)
        plt.legend()
        utils.saveplot()

        graph = GraphArtist('semilogy')
        graph.histogram(n, bins * VNS, linestyle='gray')
        self.artistplot_landau_and_gamma(graph, x, p_gamma, p_landau)
        graph.set_xlabel(r"Pulse integral [\si{\volt\nano\second}]")
        graph.set_ylabel("Count")
        graph.set_xlimits(0, 30)
        graph.set_ylimits(1e1, 1e4)
        artist.utils.save_graph(graph, dirname='plots')

    def plot_spectrum_fit_chisq(self):
        global integrals

        if 'integrals' not in globals():
            events = self.data.root.hisparc.cluster_kascade.station_601.events
            integrals = events.col('integrals')[:, 0]

        bins = np.linspace(0, RANGE_MAX, N_BINS + 1)
        n, bins = np.histogram(integrals, bins=bins)
        x = (bins[:-1] + bins[1:]) / 2

        p_gamma, p_landau = self.full_spectrum_fit(x, n, (1., 1.),
                                                   (5e3 / .32, 3.38 / 5000, 1.))
        print "FULL FIT"
        print p_gamma, p_landau

        print "charged fraction:", self.calc_charged_fraction(x, n, p_gamma, p_landau)
        landaus = scintillator.conv_landau_for_x(x, *p_landau)
        gammas = self.gamma_func(x, *p_gamma)
        fit = landaus + gammas

        x_trunc = x.compress((LOW <= x) & (x < HIGH))
        n_trunc = n.compress((LOW <= x) & (x < HIGH))
        fit_trunc = fit.compress((LOW <= x) & (x < HIGH))

        chisq, pvalue = stats.chisquare(n_trunc, fit_trunc, ddof=5)
        chisq /= (len(n_trunc) - 1 - 5)
        print "Chi-square statistic:", chisq, pvalue

        plt.figure()

        plt.plot(x * VNS, n)
        self.plot_landau_and_gamma(x, p_gamma, p_landau)
        #plt.plot(x_trunc * VNS, fit_trunc, linewidth=4)

        plt.axvline(LOW * VNS)
        plt.axvline(HIGH * VNS)

        plt.xlabel("Pulse integral [V ns]")
        plt.ylabel("Count")
        plt.yscale('log')
        plt.xlim(0, 20)
        plt.ylim(1e2, 1e5)
        plt.title(r"$\chi^2_{red}$: %.2f, p-value: %.2e" % (chisq, pvalue))
        utils.saveplot()

        plt.figure()
        plt.plot(x_trunc * VNS, n_trunc - fit_trunc)
        plt.axhline(0)
        plt.xlabel("Pulse integral [V ns]")
        plt.ylabel("Data - Fit")
        plt.title(r"$\chi^2_{red}$: %.2f, p-value: %.2e" % (chisq, pvalue))
        utils.saveplot(suffix='residuals')

    def plot_landau_and_gamma(self, x, p_gamma, p_landau):
        gammas = self.gamma_func(x, *p_gamma)
        gamma_trunc = np.where(x * VNS <= 21, gammas, 0.)

        plt.plot(x * VNS, gamma_trunc, label='gamma')

        landaus = self.scintillator.conv_landau_for_x(x, *p_landau)
        plt.plot(x * VNS, landaus, label='landau/gauss')

        plt.plot(x * VNS, gamma_trunc + landaus, label='gamma + landau/gauss')

    def artistplot_landau_and_gamma(self, graph, x, p_gamma, p_landau):
        gammas = self.gamma_func(x, *p_gamma)
        gamma_trunc = np.where(x * VNS <= 21, gammas, 1e-99)

        graph.plot(x * VNS, gamma_trunc, mark=None, linestyle='dashed')

        landaus = self.scintillator.conv_landau_for_x(x, *p_landau)
        graph.plot(x * VNS, landaus, mark=None, linestyle='dashdotted')

        graph.plot(x * VNS, gamma_trunc + landaus, mark=None)

    def artistplot_alt_landau_and_gamma(self, graph, x, p_gamma, p_landau):
        gammas = self.gamma_func(x, *p_gamma)
        gamma_trunc = np.where(x * VNS <= 21, gammas, 1e-99)

        graph.plot(x * VNS, gamma_trunc, mark=None, linestyle='dashed,gray')

        landaus = self.scintillator.conv_landau_for_x(x, *p_landau)
        graph.plot(x * VNS, landaus, mark=None, linestyle='dashdotted,gray')

    def fit_gammas_to_data(self, x, y, p0):
        condition = (LOW <= x) & (x < 2000)
        x_trunc = x.compress(condition)
        y_trunc = y.compress(condition)
        popt, pcov = optimize.curve_fit(self.gamma_func, x_trunc, y_trunc,
                                        p0=p0, sigma=np.sqrt(y_trunc))
        return popt

    def gamma_func(self, x, N, a):
        return N * x ** -a

    def fit_conv_landau_to_data(self, x, y, p0):
        popt = optimize.fmin(self.scintillator.residuals, p0,
                             (x, y, 4500, 5500), disp=0)
        return popt

    def fit_complete(self, x, y, p_gamma, p_landau):
        p0 = list(p_gamma) + list(p_landau)
        popt = optimize.fmin(self.complete_residuals, p0,
                             (self.scintillator, x, y, LOW, HIGH),
                             maxfun=100000, disp=0)
        return popt[:2], popt[2:]

    def constrained_fit_complete(self, x, y, p_gamma, p_landau):
        N_gamma = p_gamma[0]
        N_landau = p_landau[0]
        popt = optimize.fmin(self.constrained_complete_residuals,
                             (N_gamma, N_landau),
                             (self.scintillator, x, y, p_gamma,
                              p_landau, LOW, HIGH),
                             maxfun=100000, disp=0)
        p_gamma[0] = popt[0]
        p_landau[0] = popt[1]
        return p_gamma, p_landau

    def complete_residuals(self, par, scintillator, x, y, a, b):
        landaus = scintillator.conv_landau_for_x(x, *par[2:])
        gammas = self.gamma_func(x, *par[:2])
        y_exp = landaus + gammas

        y_trunc = y.compress((a <= x) & (x < b))
        y_exp_trunc = y_exp.compress((a <= x) & (x < b))

        # Make sure no zeroes end up in denominator of chi_squared
        y_trunc = np.where(y_trunc != 0., y_trunc, 1.)

        chisquared = ((y_trunc - y_exp_trunc) ** 2 / y_trunc).sum()
        return chisquared

    def constrained_complete_residuals(self, par, scintillator, x, y,
                                       p_gamma, p_landau, a, b):
        full_par = (par[0], p_gamma[1], par[1], p_landau[1], p_landau[2])
        return self.complete_residuals(full_par, scintillator, x, y, a, b)

    def get_integrals_and_densities(self):
        hisparc = self.data.root.hisparc.cluster_kascade.station_601.events
        kascade = self.data.root.kascade.events
        c_index = self.data.root.kascade.c_index
        h_index = c_index.col('h_idx')
        k_index = c_index.col('k_idx')

        intg = hisparc.read_coordinates(h_index, 'integrals')[:, 0]

        dens_e = kascade.read_coordinates(k_index, 'dens_e')[:, 0]
        dens_mu = kascade.read_coordinates(k_index, 'dens_mu')[:, 0]
        theta = kascade.read_coordinates(k_index, 'zenith')
        dens = dens_e + dens_mu
        dens_on_ground = dens * np.cos(theta)

        return intg, dens_on_ground

    def full_fit_on_data(self, integrals, p0):
        bins = np.linspace(0, RANGE_MAX, N_BINS + 1)
        n, bins = np.histogram(integrals, bins=bins)
        x = (bins[:-1] + bins[1:]) / 2

        p_gamma, p_landau = self.full_spectrum_fit(x, n, p0[:2], p0[2:])
        return list(p_gamma) + list(p_landau)

    def determine_charged_fraction(self, integrals, p0):
        bins = np.linspace(0, RANGE_MAX, N_BINS + 1)
        n, bins = np.histogram(integrals, bins=bins)
        x = (bins[:-1] + bins[1:]) / 2

        p_gamma, p_landau = self.constrained_full_spectrum_fit(x, n, p0[:2], p0[2:])
        return self.calc_charged_fraction(x, n, p_gamma, p_landau)

    def plot_detection_efficiency(self):
        integrals, dens = self.get_integrals_and_densities()

        popt = self.full_fit_on_data(integrals,
                                      (1., 1., 5e3 / .32, 3.38 / 5000, 1.))

        x, y, yerr = [], [], []
        dens_bins = np.linspace(0, 10, 51)
        for low, high in zip(dens_bins[:-1], dens_bins[1:]):
            sel = integrals.compress((low <= dens) & (dens < high))
            x.append((low + high) / 2)
            frac = self.determine_charged_fraction(sel, popt)
            y.append(frac)
            yerr.append(np.sqrt(frac * len(sel)) / len(sel))
            print (low + high) / 2, len(sel)
            self.plot_full_spectrum_fit_in_density_range(sel, popt, low, high)
        print

        plt.figure()
        plt.errorbar(x, y, yerr, fmt='o', label='data', markersize=3.)

        popt, pcov = optimize.curve_fit(self.conv_p_detection, x, y, p0=(1.,))
        print "Sigma Gauss:", popt

        x2 = plt.linspace(0, 10, 101)
        plt.plot(x2, self.p_detection(x2), label='poisson')
        plt.plot(x2, self.conv_p_detection(x2, *popt), label='poisson/gauss')

        plt.xlabel("Charged particle density [$m^{-2}$]")
        plt.ylabel("Detection probability")
        plt.ylim(0, 1.)
        plt.legend(loc='best')
        utils.saveplot()

        graph = GraphArtist()
        graph.plot(x2, self.p_detection(x2), mark=None)
        graph.plot(x2, self.conv_p_detection(x2, *popt), mark=None,
                   linestyle='dashed')
        graph.plot(x, y, yerr=yerr, linestyle=None)
        graph.set_xlabel(
            r"Charged particle density [\si{\per\square\meter}]")
        graph.set_ylabel("Detection probability")
        graph.set_xlimits(min=0)
        graph.set_ylimits(min=0)
        artist.utils.save_graph(graph, dirname='plots')

    def plot_full_spectrum_fit_in_density_range(self, sel, popt, low, high):
        bins = np.linspace(0, RANGE_MAX, N_BINS + 1)
        n, bins = np.histogram(sel, bins=bins)
        x = (bins[:-1] + bins[1:]) / 2

        p_gamma, p_landau = self.constrained_full_spectrum_fit(x, n, popt[:2], popt[2:])

        plt.figure()
        plt.plot(x * VNS, n, label='data')
        self.plot_landau_and_gamma(x, p_gamma, p_landau)

        y_charged = self.calc_charged_spectrum(x, n, p_gamma, p_landau)
        plt.plot(x * VNS, y_charged, label='charged particles')

        plt.yscale('log')
        plt.xlim(0, 50)
        plt.ylim(ymin=1)
        plt.xlabel("Pulse integral [V ns]")
        plt.ylabel("Count")
        plt.legend()
        suffix = '%.1f-%.1f' % (low, high)
        suffix = suffix.replace('.', '_')
        utils.saveplot(suffix)

        n = np.where(n > 0, n, 1e-99)
        y_charged = np.where(y_charged > 0, y_charged, 1e-99)

        graph = GraphArtist('semilogy')
        graph.histogram(n, bins * VNS, linestyle='gray')
        self.artistplot_alt_landau_and_gamma(graph, x, p_gamma, p_landau)
        graph.histogram(y_charged, bins * VNS)
        graph.set_xlabel(r"Pulse integral [\si{\volt\nano\second}]")
        graph.set_ylabel("Count")
        graph.set_title(r"$\SI{%.1f}{\per\square\meter} \leq \rho_\mathrm{charged}$ < $\SI{%.1f}{\per\square\meter}$" % (low, high))
        graph.set_xlimits(0, 30)
        graph.set_ylimits(1e0, 1e4)
        artist.utils.save_graph(graph, suffix, dirname='plots')

    p_detection = np.vectorize(lambda x: 1 - np.exp(-.5 * x) if x >= 0 else 0.)

    def conv_p_detection(self, x, sigma):
        x_step = x[-1] - x[-2]
        x2 = np.arange(-2 * max(x), 2 * max(x) + x_step / 2, x_step)
        g = stats.norm(scale=sigma).pdf
        y2 = landau.discrete_convolution(self.p_detection, g, x2)
        y = np.interp(x, x2, y2)
        return y


if __name__ == '__main__':
    np.seterr(invalid='ignore', divide='ignore')

    if 'data' not in globals():
        data = tables.open_file('kascade.h5', 'r')

    utils.set_prefix('EFF-')
    artist.utils.set_prefix('EFF-')
    efficiency = ReconstructionEfficiency(data)
    efficiency.main()
