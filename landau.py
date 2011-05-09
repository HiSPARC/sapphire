""" Landau distribution function

    This module computes the Landau distribution, which governs the
    fluctuations in energy loss of particles travelling through a
    relatively thin layer of matter.

    Currently, this module only contains functions to calculate the exact
    function using two integral representations of the defining complex
    integral.  This should be extended by approximations when the need for
    doing serious work arises.

"""
import numpy as  np
from numpy import pi, Inf, sin, cos, exp, log, sqrt, arctan, vectorize, \
                  convolve
from scipy import integrate, stats

#def residuals(p, y, x):
#    N, location, scale = p
#    return y - N * pdf(x, location, scale)

@vectorize
def pdf(lf):
    if lf < -10:
        return 0.
    elif lf < 0:
        sf = exp(-lf - 1)
        integrant = integrate.quad(pdf_kernel, 0, Inf, args=(sf,))[0]
        return 1 / pi * exp(-sf) * integrant
    else:
        integrant = integrate.quad(pdf_kernel2, 0, Inf, args=(lf,))[0]
        return 1 / pi * integrant

def pdf_kernel(y, sf):
    return (exp(sf / 2 * log(1 + y ** 2 / sf ** 2) - y * arctan(y / sf)) *
            cos(.5 * y * log(1 + y ** 2 / sf ** 2) - y + sf * arctan(y / sf)))

def pdf_kernel2(u, lf):
    return exp(-lf * u) * u ** -u * sin(pi * u)


class Scintillator:
    thickness = .02 # m
    xi        = 0.172018 # MeV
    epsilon   = 3.10756e-11
    delta     = 2.97663
    Euler     = 0.577215665

    _lf0 = log(xi) - log(epsilon) + 1 - Euler - delta

    def landau_pdf(self, Delta):
        lf = Delta / self.xi - self._lf0
        return pdf(lf) / self.xi

    def conv_landau(self, x, mev_scale=1, count_scale=1, gauss_scale=1):
        f = self.landau_pdf
        g = stats.norm(scale=gauss_scale).pdf
        return count_scale * discrete_convolution(f, g, mev_scale * x)

    def residuals(self, p, xdata, ydata, a, b):
        mev_scale, count_scale, gauss_scale = p
        yfit = self.conv_landau(xdata, mev_scale, count_scale,
                                gauss_scale)

        yfit = yfit.compress((a <= xdata) & (xdata < b))
        ydata = ydata.compress((a <= xdata) & (xdata < b))
        weights = 1. / ydata
        return ((weights * (yfit - ydata)) ** 2).sum()
        #return ((yfit.compress((a <= xdata) & (xdata < b)) -
        #        ydata.compress((a <= xdata) & (xdata < b))) ** 2).sum()


@vectorize
def fixed_convolution(f, g, t, a=-10, b=10, N=50):
    func = lambda tau: f(tau) * g(t - tau)
    tau = linspace(a, b, N)
    return integrate.trapz(func(tau), tau)

@vectorize
def convolution(f, g, t):
    func = lambda tau: f(tau) * g(t - tau)
    return integrate.quad(func, -integrate.Inf, integrate.Inf)[0]

def discrete_convolution(f, g, t):
    if -min(t) != max(t):
        raise RuntimeError("Range needs to be symmetrical around zero.")

    dt = t[1] - t[0]
    return dt * convolve(f(t), g(t), mode='same')
