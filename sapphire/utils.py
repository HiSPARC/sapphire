"""Utilities

The module contains some commonly functions and classes.

"""
import numpy as np
from progressbar import ProgressBar, ETA, Bar, Percentage


def pbar(iterable, length=None, show=True, **kwargs):
    """Get a new progressbar with our default widgets

    :param iterable: the iterable over which will be looped.
    :param length: in case iterable is a generator, this should be its
                   expected length.
    :param show: boolean, if False simply return the iterable.
    :returns: a new iterable which iterates over the same elements as
              the input, but shows a progressbar if possible.

    """
    if not show:
        return iterable

    if length is None:
        try:
            length = len(iterable)
        except TypeError:
            pass

    if length:
        pb = ProgressBar(widgets=[ETA(), Bar(), Percentage()], maxval=length,
                         **kwargs)
        return pb(iterable)
    else:
        return iterable


def ceil_in_base(value, base):
    """Get nearest multiple of base above the value"""

    return base * np.ceil(value / base)


def floor_in_base(value, base):
    """Get nearest multiple of base below the value"""

    return base * np.floor(value / base)


def gauss(x, N, mu, sigma):
    """Gaussian distribution

    To be used for fitting where the integral is not 1.

    """
    return N * norm.pdf(x, mu, sigma)
