""" Perform various coordinate transformations

    This module performs various coordinate transformations, based on some
    well-known formulas.

"""
from math import sin, cos, sqrt, radians

from numpy import matrix


class WGS84Datum(object):
    """Definition of the WGS84 datum

    These definitions are taken from
    http://en.wikipedia.org/wiki/Geodetic_system, believing that enough
    editors have gone over them to make sure they are correct.

    """
    # Defining constants
    a = 6378137.
    f = 1 / 298.257223563

    # Derived constants
    b = a * (1 - f)
    e = sqrt(2 * f - f ** 2)
    eprime = sqrt(f * (2 - f) / (1 - f) ** 2)


class FromWGS84ToENUTransformation(object):
    geode = WGS84Datum()

    def __init__(self, ref_llacoordinates):
        self.ref_lla = ref_llacoordinates
        self.ref_XYZ = self.lla_to_ecef(ref_llacoordinates)

    def transform(self, coordinates):
        """Transfrom WGS84 coordinates to ENU coordinates"""

        return self.ecef_to_enu(self.lla_to_ecef(coordinates))

    def lla_to_ecef(self, coordinates):
        """Convert from LLA coordinates to ECEF coordinates

        LLA:  Latitude, Longitude, Altitude
        ECEF: Earth-Centered, Earth-Fixed

        The conversion formulas are taken from
        http://en.wikipedia.org/wiki/Geodetic_system#From_geodetic_to_ECEF
        but slightly reworked.

        Mind that the input is expected to be in degrees, as is standard in
        coordinate notation.

        :param coordinates: tuple of latitude, longitude (both in degrees)
                            and altitude (in meters).

        """
        latitude, longitude, altitude = coordinates

        latitude = radians(latitude)
        longitude = radians(longitude)

        a = self.geode.a
        b = self.geode.b
        e = self.geode.e

        N = a / sqrt(1 - e ** 2 * sin(latitude) ** 2)

        X = (N + altitude) * cos(latitude) * cos(longitude)
        Y = (N + altitude) * cos(latitude) * sin(longitude)
        Z = (b ** 2 / a ** 2 * N + altitude) * sin(latitude)

        return X, Y, Z

    def ecef_to_enu(self, coordinates):
        """Convert from ECEF coordinates to ENU coordinates

        ECEF: Earth-Centered, Earth-Fixed
        ENU: East, North, Up

        The conversion formulas are taken from
        http://en.wikipedia.org/wiki/Geodetic_system#From_ECEF_to_ENU

        :param coordinates: a tuple containing the ECEF coordinates (in meters)
                            of the point to transform

        """
        latitude, longitude, altitude = self.ref_lla
        Xr, Yr, Zr = self.ref_XYZ
        X, Y, Z = coordinates

        lat = radians(latitude)
        lon = radians(longitude)

        transformation = matrix([
            [           -sin(lon),             cos(lon),       0.],
            [-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat)],
            [ cos(lat) * cos(lon),  cos(lat) * sin(lon), sin(lat)]])

        coordinates = matrix([[X - Xr], [Y - Yr], [Z - Zr]])

        return (transformation * coordinates).A1
