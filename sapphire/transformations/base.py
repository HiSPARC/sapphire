"""Perform base conversions

Currently supports conversion between base 10 (decimal) and
base 60 (sexagesimal).

"""
from numpy import modf


def decimal_to_sexagesimal(decimal):
    """Convert decimal hours or degrees to sexagesimal.

    :param decimal: decimal number to be converted to sexagismal.
    :returns: tuple of either (hours, minutes, seconds) or
              (degrees, arcminutes, arcseconds)

    """
    fractional, integral = modf(decimal)
    min_fractional, minutes = modf(fractional * 60)
    seconds = min_fractional * 60.
    return integral.astype(int), minutes.astype(int), seconds


def sexagesimal_to_decimal(hd, minutes, seconds):
    """Convert sexagesimal hours or degrees to decimal.

    :param hd: hours or degrees.
    :param minutes: minutes or arcminutes.
    :param seconds: seconds or arcseconds.
    :returns: decimal hours or degrees.

    """
    return hd + minutes / 60. + seconds / 3600.
