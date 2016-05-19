"""Utilities

The module contains some commonly functions and classes.

"""
from __future__ import division

from functools import wraps
from bisect import bisect_right
from distutils.spawn import find_executable

from numpy import floor, ceil, round, arcsin, sin, pi, sqrt
from scipy.stats import norm
from progressbar import ProgressBar, ETA, Bar, Percentage


#: Error values used to indicate missing or bad data.
#: Code -999 is used if the reconstruction of a quantity failed.
#: Code -1 is used if that detector/sensor is not present.
ERR = [-1, -999]

#: Speed of light in vacuum in m / ns.
c = 0.299792458


def pbar(iterable, length=None, show=True, **kwargs):
    """Get a new progressbar with our default widgets

    :param iterable: the iterable over which will be looped.
    :param length: in case iterable is a generator, this should be its
                   expected length.
    :param show: boolean, if False simply return the iterable.
    :return: a new iterable which iterates over the same elements as
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
        pb = ProgressBar(widgets=[Percentage(), Bar(), ETA()], maxval=length,
                         **kwargs)
        return pb(iterable)
    else:
        return iterable


def ceil_in_base(value, base):
    """Get nearest multiple of base above the value"""

    return base * ceil(value / base)


def floor_in_base(value, base):
    """Get nearest multiple of base below the value"""

    return base * floor(value / base)


def round_in_base(value, base):
    """Get nearest multiple of base to the value"""

    return base * round(value / base)


def closest_in_list(value, items):
    """Get nearest item from a list of items to the value"""

    return min(items, key=lambda x: abs(x - value))


def get_active_index(values, value):
    """Get the index where the value fits.

    :param values: sorted list of values (e.g. list of timestamps).
    :param value: value for which to find the position (e.g. a timestamp).
    :return: index into the values list.

    """
    idx = bisect_right(values, value, lo=0)
    if idx == 0:
        idx = 1
    return idx - 1


def gauss(x, N, mu, sigma):
    """Gaussian distribution

    To be used for fitting where the integral is not 1.

    """
    return N * norm.pdf(x, mu, sigma)


def norm_angle(angle):
    """Normalize an angle to the range [-pi, pi)

    We use the range from -pi upto but not including pi to represent
    angles.

    """
    return (angle + pi) % (2 * pi) - pi


def angle_between(zenith1, azimuth1, zenith2, azimuth2):
    """Calculate the angle between two (zenith, azimuth) coordinates

    Using the haversine formula,
    from: http://www.movable-type.co.uk/scripts/latlong.html

    :param zenith#: Zenith parts of the coordinates, in radians (0, pi/2).
    :param azimuth#: Azimuth parts of the coordinates, in radians (-pi, pi).
    :return: Angle between the two coordinates.

    """
    dlat = zenith1 - zenith2
    dlon = azimuth2 - azimuth1
    a = (sin(dlat / 2) ** 2 + sin(zenith1) * sin(zenith2) * sin(dlon / 2) ** 2)
    angle = 2 * arcsin(sqrt(a))

    return angle


def distance_between(x1, y1, x2, y2):
    """Calculate the distance between two (x, y) coordinates

    :param x#: x parts of the coordinates.
    :param y#: y parts of the coordinates.
    :return: distance between the two coordinates.

    """
    d = sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    return d


def make_relative(x):
    """Make first element the origin and make rest relative to it."""

    return [xi - x[0] for xi in x]


def which(program):
    """Check if a command line program is available

    An Exception is raised if the program is not available.

    :param program: name or program to check for, e.g. 'wget'.

    """
    path = find_executable(program)
    if not path:
        raise Exception('The program %s is not available.' % program)


def memoize(obj):
    """ Memoisation cache decorator

    """
    cache = obj.cache = {}

    @wraps(obj)
    def memoizer(*args, **kwargs):
        key = str(args) + str(kwargs)
        if key not in cache:
            cache[key] = obj(*args, **kwargs)
        return cache[key]
    return memoizer
