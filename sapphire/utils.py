"""Utilities

The module contains some commonly functions and classes.

"""
from numpy import floor, ceil, arccos, cos, sin, pi
from scipy.stats import norm
from progressbar import ProgressBar, ETA, Bar, Percentage


# Error values used to indicate missing or bad data.
ERR = [-1, -999]


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

    return base * ceil(value / base)


def floor_in_base(value, base):
    """Get nearest multiple of base below the value"""

    return base * floor(value / base)


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

    Using the spherical law of cosines,
    from: http://www.movable-type.co.uk/scripts/latlong.html#cosine-law

    :param zenith#: Zenith parts of the coordinates, in radians (0, pi/2).
    :param azimuth#: Azimuth parts of the coordinates, in radians (-pi, pi).
    :returns: Angle between the two coordinates.

    """
    lat1 = pi / 2 - zenith1
    lat2 = pi / 2 - zenith2
    return arccos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) *
                  cos(azimuth1 - azimuth2))
