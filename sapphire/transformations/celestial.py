""" Perform various Celestial coordinate transformations

    This module performs transformations between different
    Celestial coordinate systems.

    Legacy transformations (all those not marked astropy):
    Formulae from: Duffett-Smith1990
    'Astronomy with your personal computer'
    ISBN 0-521-38995-X

    New transformations have been added with _astropy added to function name
    They are very exact, in the order of arcsec.
    Ethan van Woerkom is the author of the new transformations; contact him
    for further information.
"""
import datetime
import warnings

import numpy as np

from numpy import arccos, arcsin, around, cos, pi, radians, sin

from ..utils import norm_angle
from . import angles, clock


def zenithazimuth_to_equatorial(latitude, longitude, timestamp, zenith,
                                azimuth):
    """Convert Zenith Azimuth to Equatorial coordinates (J2000.0)

    :param latitude,longitude: Position of the observer on Earth in degrees.
                               North and east positive.
    :param timestamp: GPS timestamp of the observation in seconds.
    :param zenith: zenith is the angle relative to the Zenith in radians.
    :param azimuth: azimuth angle of the observation in radians.

    :return: Right ascension (ra) and Declination (dec) in radians.

    From Duffett-Smith1990, 1500 EQHOR and 1600 HRANG

    """
    altitude, alt_azimuth = zenithazimuth_to_horizontal(zenith, azimuth)
    lst = clock.gps_to_lst(timestamp, longitude)
    ra, dec = horizontal_to_equatorial(latitude, lst, altitude, alt_azimuth)

    return ra, dec


def zenithazimuth_to_horizontal(zenith, azimuth):
    """Convert from Zenith Azimuth to Horizontal coordinates

    :param zenith: Zenith in radians
    :param azimuth: Azimuth in radians
    :return altitude, alt_azimuth: Alt, Az in radians

    Zenith Azimuth is the coordinate system used by HiSPARC. Zenith is
    the angle between the zenith and the direction. Azimuth is the angle
    in the horizontal plane, from East to North (ENWS).

    Horizontal is the coordinate system as described in
    Duffett-Smith1990 p38. Altitude is the angle above the horizon and
    Azimuth the angle in the horizontal plane, from North to East (NESW).

    """
    altitude = norm_angle(pi / 2. - zenith)
    alt_azimuth = norm_angle(pi / 2. - azimuth)

    return altitude, alt_azimuth


def horizontal_to_zenithazimuth(altitude, alt_azimuth):
    """Inverse of zenithazimuth_to_horizontal is the same transformation"""

    return zenithazimuth_to_horizontal(altitude, alt_azimuth)


def horizontal_to_equatorial(latitude, lst, altitude, alt_azimuth):
    """Convert Horizontal to Equatorial coordinates (J2000.0)

    :param latitude: Position of the observer on Earth in degrees.
                     North positive.
    :param lst: Local Siderial Time observer at the time of observation
                in decimal hours.
    :param altitude: altitude is the angle above the horizon in radians.
    :param alt_azimuth: Azimuth angle in horizontal plane in radians.

    :return: Right ascension (ra) and Declination (dec) in radians.

    Warning: Inexact transformation; astropy functions preferred.

    From Duffett-Smith1990, 1500 EQHOR and 1600 HRANG

    """
    ha, dec = horizontal_to_hadec(latitude, altitude, alt_azimuth)
    ra = ha_to_ra(ha, lst)

    return ra, dec


def horizontal_to_hadec(latitude, altitude, alt_azimuth):
    """Convert Horizontal to Hour Angle and Declination

    :param latitude: Position of the observer on Earth in degrees.
                     North positive.
    :param altitude: altitude is the angle above the horizon in radians.
    :param alt_azimuth: Azimuth angle in horizontal plane in radians.

    :return: Hour angle (ha) and Declination (dec) in radians.

    Warning: Inexact transformation; astropy functions preferred.

    From Duffett-Smith1990, 1500 EQHOR and 1600 HRANG

    """

    slat = sin(radians(latitude))
    clat = cos(radians(latitude))
    sazi = sin(alt_azimuth)
    cazi = cos(alt_azimuth)
    salt = sin(altitude)
    calt = cos(altitude)

    dec = arcsin((salt * slat) + (calt * clat * cazi))

    # Round to prevent value beyond allowed range for arccos.
    cha = around((salt - (slat * sin(dec))) / (clat * cos(dec)), 15)
    ha = arccos(cha)

    if sazi > 0:
        ha = 2 * pi - ha

    return ha, dec


def ha_to_ra(ha, lst):
    """Convert Hour angle to right ascension

    :param ha: Hour angle in radians.
    :param lst: Local Siderial Time observer at the time of observation
                in decimal hours.

    :return: Right ascension (ra) in radians.

    """
    ra = (angles.hours_to_radians(lst) - ha)
    ra %= 2 * pi

    return ra


def equatorial_to_zenithazimuth(latitude, longitude, timestamp,
                                right_ascension, declination):
    """Convert Equatorial (J2000.0) to Zenith Azimuth coordinates

    :param latitude,longitude: Position of the observer on Earth in degrees.
                               North and east positive.
    :param timestamp: GPS timestamp of the observation in seconds.
    :param right_ascension: right_ascension of the observation in radians.
    :param declination: declination of the observation in radians.

    :return: zenith and azimuth in radians.

    This function was renamed from equatorial_to_horizontal to
    equatorial_to_zenithazimuth in order to make it operate as the name does.

    From Duffett-Smith1990, 1500 EQHOR and 1600 HRANG

    """
    lst = clock.gps_to_lst(timestamp, longitude)
    ha = (angles.hours_to_radians(lst) - right_ascension)
    ha %= 2 * pi

    slat = sin(radians(latitude))
    clat = cos(radians(latitude))
    sha = sin(ha)
    cha = cos(ha)
    sdec = sin(declination)
    cdec = cos(declination)

    altitude = arcsin((sdec * slat) + (cdec * clat * cha))
    alt_azimuth = arccos((sdec - (slat * sin(altitude))) /
                         (clat * cos(altitude)))

    if sha > 0:
        alt_azimuth = 2 * pi - alt_azimuth

    zenith, azimuth = horizontal_to_zenithazimuth(altitude, alt_azimuth)

    return zenith, azimuth


try:
    # This try-except structure has been implemented,
    # to accommodate those without astropy.
    import astropy.units as u

    from astropy.coordinates import EarthLocation, SkyCoord
    from astropy.time import Time

    def zenithazimuth_to_equatorial_astropy(latitude, longitude, utc_timestamp,
                                            zenaz_coordinates):
        """ Converts iterables of tuples of zenithazimuth
            to equatorial coordinates

        :param latitude: Latitude in decimal degrees
        :param longitude: Longitude in decimal degrees
        :param utc_timestamp: Unix UTC timestamp integer
        :param zenaz_coordinates: np.array of tuples (zen, az) in radians
        :return: np.array of tuples (ra, dec) in radians

        For increased speed using array input is recommended.
        """

        # Convert and flip order of zenaz coordinates, done in numpy for speed
        zenaz_coordinates = np.array(zenaz_coordinates)
        zenaz_coordinates = 0.5 * np.pi - zenaz_coordinates
        horizontal_coordinates = np.unwrap(zenaz_coordinates[:, [1, 0]])

        # Normalise angle
        horizontal_coordinates = norm_angle(horizontal_coordinates)

        return horizontal_to_equatorial_astropy(latitude, longitude,
                                                utc_timestamp,
                                                horizontal_coordinates)

    def equatorial_to_zenithazimuth_astropy(latitude, longitude,
                                            utc_timestamp,
                                            equatorial_coordinates):
        """ Converts iterables of tuples of equatorial
            to zenithazimuth coordinates

        :param latitude: Latitude in decimal degrees
        :param longitude: Longitude in decimal degrees
        :param utc_timestamp: Unix UTC timestamp integer
        :param equatorial_coordinates: np.array of tuples (ra, dec) in radians
        :return: np.array of tuples (zen, az) in radians

        For increased speed using array input is recommended.
        """

        equatorial_coordinates = np.array(equatorial_coordinates)
        horizontal_coordinates = equatorial_to_horizontal_astropy(
            latitude, longitude, utc_timestamp, equatorial_coordinates)

        # Convert and flip order of zenaz coordinates, done in numpy for speed
        horizontal_coordinates = np.array(horizontal_coordinates)
        horizontal_coordinates = 0.5 * np.pi - horizontal_coordinates
        zenaz_coordinates = horizontal_coordinates[:, [1, 0]]

        # Normalise angle
        zenaz_coordinates = norm_angle(zenaz_coordinates)

        return zenaz_coordinates

    def equatorial_to_horizontal_astropy(latitude, longitude,
                                         utc_timestamp,
                                         equatorial_coordinates):
        """ Converts iterables of tuples of equatorial coordinates
            to horizontal coordinates

        :param latitude: Latitude in decimal degrees
        :param longitude: Longitude in decimal degrees
        :param utc_timestamp: Unix UTC timestamp integer
        :param equatorial_coordinates: np.array of tuples (ra, dec) in radians
        :return: np.array of tuples (az, alt) in radians

        For increased speed using array input is recommended.
        """
        # For speed in numpy
        equatorial_coordinates = np.array(equatorial_coordinates)

        location = EarthLocation(longitude, latitude)
        t = Time(datetime.datetime.utcfromtimestamp(utc_timestamp))
        equatorial_frame = SkyCoord(equatorial_coordinates, location=location,
                                    obstime=t, unit=u.rad, frame='icrs')
        horizontal_frame = equatorial_frame.transform_to('altaz')

        return np.array((horizontal_frame.az.rad, horizontal_frame.alt.rad)).T

    def horizontal_to_equatorial_astropy(latitude, longitude,
                                         utc_timestamp,
                                         horizontal_coordinates):
        """ Converts iterables of tuples of
            horizontal coordinates to equatorial coordinates

        :param latitude: Latitude in decimal degrees
        :param longitude: Longitude in decimal degrees
        :param utc_timestamp: Unix UTC timestamp integer
        :param horizontal_coordinates: np.array of tuples (az, alt) in radians
        :return: np.array of tuples (ra, dec) in radians
        """
        # For speed in numpy
        horizontal_coordinates = np.array(horizontal_coordinates)

        location = EarthLocation(longitude, latitude)
        t = Time(datetime.datetime.utcfromtimestamp(utc_timestamp))
        horizontal_frame = SkyCoord(horizontal_coordinates, location=location,
                                    obstime=t, unit=u.rad, frame='altaz')
        equatorial_frame = horizontal_frame.transform_to('icrs')

        return np.array((equatorial_frame.ra.rad, equatorial_frame.dec.rad)).T


except ImportError as e:
    warnings.warn(str(e) + "\nImport of astropy failed", ImportWarning)
