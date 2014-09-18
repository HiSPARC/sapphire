""" Perform various Celestial coordinate transformations

    This module performs transformations between different
    Celestial coordinate systems.

    Formulae from: Duffett-Smith1990
    'Astronomy with your personal computer'
    ISBN 0-521-38995-X

"""
from numpy import (arcsin, arccos, arctan2, cos, sin,
                   array, radians, degrees, pi)

from . import clock, angles


def horizontal_to_equitorial(longitude, latitude, timestamp, azimuth, zenith):
    """Convert Horizontal to Equatorial coordinates (J2000.0)

    :param longitude,latitude: Position of the observer on Earth in degrees.
                               North and east positive.
    :param timestamp: GPS timestamp of the observation.
    :param azimuth: zenith angle of the observation in radians.
    :param zenith: azimuth angle of the observation in radians.

    :returns: Right ascension (ra) and Declination (dec) in radians.

    From Duffett-Smith1990, 1500 EQHOR and 1600 HRANG

    """
    # altitude is the angle above the horizon
    altitude = pi / 2. - zenith

    slat = sin(radians(latitude))
    clat = cos(radians(latitude))
    sazi = sin(azimuth)
    cazi = cos(azimuth)
    salt = sin(altitude)
    calt = cos(altitude)

    dec = arcsin((salt * slat) + (calt * clat * cazi))
    HA = arccos((salt - (slat * sin(dec))) / (clat * cos(dec)))

    if sazi > 0:
        HA = 2 * pi - HA

    lst = clock.gps_to_lst(timestamp, longitude)
    ra = (angles.hours_to_radians(lst) - HA)
    ra %= 2 * pi

    return ra, dec


def equitorial_to_horizontal(longitude, latitude, timestamp, right_ascension,
                             declination):
    """Convert Equatorial (J2000.0) to Horizontal coordinates

    :param longitude,latitude: Position of the observer on Earth in degrees.
                               North and east positive.
    :param timestamp: GPS timestamp of the observation.
    :param right_ascension: right_ascension of the observation in radians.
    :param declination: declination of the observation in radians.

    :returns: azimuth and zenith in radians.

    From Duffett-Smith1990, 1500 EQHOR and 1600 HRANG

    """
    lst = clock.gps_to_lst(timestamp, longitude)
    HA = (angles.hours_to_radians(lst) - right_ascension)
    HA %= 2 * pi

    slat = sin(radians(latitude))
    clat = cos(radians(latitude))
    sha = sin(HA)
    cha = cos(HA)
    sdec = sin(declination)
    cdec = cos(declination)

    altitude = arcsin((sdec * slat) + (cdec * clat * cha))
    azimuth = arccos((sdec - (slat * sin(altitude))) / (clat * cos(altitude)))

    if sha > 0:
        azimuth = 2 * pi - azimuth

    # altitude is the angle above the horizon
    zenith = pi / 2. - altitude

    return azimuth, zenith
