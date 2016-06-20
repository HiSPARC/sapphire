""" Landau distribution function

    This module computes the Landau distribution, which governs the
    fluctuations in energy loss of particles travelling through a
    relatively thin layer of matter.

    Currently, this module only contains functions to calculate the exact
    function using two integral representations of the defining complex
    integral.  This should be extended by approximations when the need for
    doing serious work arises.

    References are made to Fokkema2012, DOI: 10.3990/1.9789036534383.

"""
import warnings

from numpy import (pi, Inf, sin, cos, exp, log, arctan, vectorize,
                   convolve, linspace, interp)
from scipy import integrate, stats


@vectorize
def pdf(lf):
    """The Landau probability density function

    Fokkema2012, eq 2.13.

    """
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
    """The Landau kernel

    Fokkema2012, eq 2.13.

    """
    return exp(-lf * u) * u ** -u * sin(pi * u)


class Scintillator(object):
    thickness = .02  # m
    xi = 0.172018  # MeV, Fokkema2012, eq 2.12.
    epsilon = 3.10756e-11  # Fokkema2012, eq 2.11.
    delta = 2.97663  # Delta
    Euler = 0.577215665  # Euler-Mascheroni constant

    mev_scale = 1
    gauss_scale = 1

    # Fokkema2012, eq 2.10.
    _lf0 = log(xi) - log(epsilon) + 1 - Euler - delta

    full_domain = linspace(-100, 100, 1000)
    pdf_values = None
    pdf_domain = full_domain.compress(full_domain >= -5)

    def landau_pdf(self, Delta):
        """The Landau energy loss distribution function

        Fokkema2012, eq 2.9, where lf is eq 2.10.

        :param Delta: Energy loss in the scintillator.
        :return: energy loss probability.

        """
        lf = self.lf(Delta)
        return self.pdf(lf) / self.xi

    def lf(self, Delta):
        """Calculate the lambda parameter

        Fokkema2012, eq 2.10.
        With additional shift by delta.

        :param Delta: Energy loss in the scintillator.
        :return: lambda parameter.

        """
        return Delta / self.xi - self._lf0

    def pdf(self, lf):
        """The Landau probability density function

        Fokkema2012, eq 2.13.

        :param lf: lambda parameter.
        :return: probability.

        """
        if self.pdf_values is not None:
            return interp(lf, self.pdf_domain, self.pdf_values)
        else:
            # Generate pre-computed values for Landau PDF
            self.pdf_values = pdf(self.pdf_domain)
            return self.pdf(lf)

    def conv_landau_for_x(self, x, count_scale=1, mev_scale=None,
                          gauss_scale=None):
        """Landau convolved with Gaussian

        Fokkema2012, eq 5.4.

        :param x: energy loss(es) for which to get the probability.
        :param count_scale: total number of counts.
        :param mev_scale: number of MeV per unit of x.
        :param gauss_scale: width of the normal distribution.
        :return: probability.

        """
        if mev_scale is None:
            mev_scale = self.mev_scale
        if gauss_scale is None:
            gauss_scale = self.gauss_scale

        f = self.landau_pdf
        g = stats.norm(scale=gauss_scale).pdf
        x_domain = self.full_domain

        y_calc = count_scale * discrete_convolution(f, g, x_domain)
        x_calc = x_domain / mev_scale

        y = interp(x, x_calc, y_calc)
        return y

    def conv_landau(self, x, count_scale=1, mev_scale=None,
                    gauss_scale=None):
        """Bare-bones convoluted landau function

        This thing is fragile.  Use with great care!  First and foremost,
        x must be symmetrical around zero.  Second, x must contain most of
        the Landau function (including a significant part of the tail).
        If not, the results cannot be trusted!

        Better use conv_landau_for_x, which better handles this.

        """
        warnings.warn("Better be sure you know that you're doing!")

        if not mev_scale:
            mev_scale = self.mev_scale
        if gauss_scale is None:
            gauss_scale = self.gauss_scale

        f = self.landau_pdf
        g = stats.norm(scale=gauss_scale).pdf
        return count_scale * discrete_convolution(f, g, mev_scale * x)

    def residuals(self, p, xdata, ydata, a, b):
        count_scale, mev_scale, gauss_scale = p
        self.mev_scale = mev_scale
        self.gauss_scale = gauss_scale

        return self._residuals(xdata, ydata, mev_scale, count_scale,
                               gauss_scale, a, b)

    def constrained_residuals(self, p, xdata, ydata, a, b):
        count_scale = p
        mev_scale = self.mev_scale
        gauss_scale = self.gauss_scale

        return self._residuals(xdata, ydata, mev_scale, count_scale,
                               gauss_scale, a, b)

    def _residuals(self, xdata, ydata, mev_scale, count_scale,
                   gauss_scale, a, b):
        yfit = self.conv_landau_for_x(xdata, count_scale, mev_scale,
                                      gauss_scale)

        yfit = yfit.compress((a <= xdata) & (xdata < b))
        ydata = ydata.compress((a <= xdata) & (xdata < b))
        return ((yfit - ydata) ** 2 / ydata).sum()


def discrete_convolution(f, g, t):
    """Discrete convolution

    :param f,g: two functions that take one argument (t).
    :param t: values for which the functions will be evaluated, and the along
              which the convolution will be performed.
    :return: convolution of the two functions.

    """
    if abs(min(t) + max(t)) > 1e-6:
        raise RuntimeError("Range needs to be symmetrical around zero.")

    dt = t[1] - t[0]
    return dt * convolve(f(t), g(t), mode='same')
