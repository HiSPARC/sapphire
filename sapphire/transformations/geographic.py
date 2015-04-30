""" Perform various coordinate transformations

    This module performs various coordinate transformations, based on some
    well-known formulas.

"""
from math import sin, cos, atan2, sqrt, radians, degrees

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

        return self.lla_to_enu(coordinates)

    def lla_to_enu(self, coordinates):
        """Transfrom WGS84 coordinates to ENU coordinates"""

        return self.ecef_to_enu(self.lla_to_ecef(coordinates))

    def enu_to_lla(self, coordinates):
        """Transfrom WGS84 coordinates to ENU coordinates"""

        return self.ecef_to_lla(self.enu_to_ecef(coordinates))

    def lla_to_ecef(self, coordinates):
        """Convert from LLA coordinates to ECEF coordinates

        LLA: Latitude, Longitude, Altitude
        ECEF: Earth-Centered, Earth-Fixed

        The conversion formulas are taken from
        http://en.wikipedia.org/wiki/Geodetic_system#From_geodetic_to_ECEF
        but slightly reworked.

        Mind that the input is expected to be in degrees, as is standard in
        coordinate notation.

        :param coordinates: tuple of latitude, longitude (both in degrees)
                            and altitude (in meters).
        :return: ECEF coordinates (in meters).

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

    def ecef_to_lla(self, coordinates):
        """Convert from ECEF coordinates to LLA coordinates

        ECEF: Earth-Centered, Earth-Fixed
        LLA: Latitude, Longitude, Altitude

        The conversion formulas are taken from
        https://gist.github.com/klucar/1536054

        :param coordinates: tuple of X, Y, and Z (in meters).
        :return: latitude, longitude (in degrees) and altitude (in meters).

        """
        X, Y, Z = coordinates

        a = self.geode.a
        b = self.geode.b
        e = self.geode.e
        eprime = self.geode.eprime

        p = sqrt(X ** 2 + Y ** 2)
        th = atan2(a * Z, b * p)

        longitude = atan2(Y, X)
        latitude = atan2((Z + eprime ** 2 * b * sin(th) ** 3),
                         (p - e ** 2 * a * cos(th) ** 3))
        N = a / sqrt(1 - e ** 2 * sin(latitude) ** 2)
        altitude = p / cos(latitude) - N

        return degrees(latitude), degrees(longitude), altitude

    def ecef_to_enu(self, coordinates):
        """Convert from ECEF coordinates to ENU coordinates

        ECEF: Earth-Centered, Earth-Fixed
        ENU: East, North, Up

        The conversion formulas are taken from
        http://en.wikipedia.org/wiki/Geodetic_system#From_ECEF_to_ENU

        :param coordinates: a tuple containing the ECEF coordinates (in meters)
                            of the point to transform
        :return: east, north, and up (in meters).

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

    def enu_to_ecef(self, coordinates):
        """Convert from ENU coordinates to ECEF coordinates

        ENU: East, North, Up
        ECEF: Earth-Centered, Earth-Fixed

        :param coordinates: a tuple containing the ENU coordinates (in meters).
        :return: ECEF coordinates (in meters).

        """
        latitude, longitude, altitude = self.ref_lla
        Xr, Yr, Zr = self.ref_XYZ

        lat = radians(latitude)
        lon = radians(longitude)

        transformation = matrix([
            [-sin(lon), -sin(lat) * cos(lon), cos(lat) * cos(lon)],
            [ cos(lon), -sin(lat) * sin(lon), cos(lat) * sin(lon)],
            [       0.,             cos(lat),            sin(lat)]])

        x, y, z = (transformation * matrix(coordinates).T).A1

        return x + Xr, y + Yr, z + Zr
