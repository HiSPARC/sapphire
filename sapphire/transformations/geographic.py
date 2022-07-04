""" Perform various coordinate transformations

    This module performs various coordinate transformations, based on some
    well-known formulas.

"""
from math import atan2, cos, degrees, radians, sin, sqrt

from numpy import array


class WGS84Datum:
    """Definition of the WGS84 datum

    These definitions are taken from
    https://en.wikipedia.org/wiki/Geodetic_system, believing that enough
    editors have gone over them to make sure they are correct.

    """
    # Defining constants
    a = 6378137.
    f = 1 / 298.257223563

    # Derived constants
    b = a * (1 - f)
    e = sqrt(2 * f - f ** 2)
    eprime = sqrt(f * (2 - f) / (1 - f) ** 2)


class FromWGS84ToENUTransformation:

    """Convert between various geographic coordinate systems

    This class converts coordinates between LLA, ENU, and ECEF.

    """

    geode = WGS84Datum()

    def __init__(self, ref_llacoordinates):
        """Initialize the transformation object.

        :param ref_llacoordinates: reference latitude, longitude, and altitude
            coordinates. These are used as origin for ENU coordinates.

        """
        self.ref_lla = ref_llacoordinates
        self.ref_ecef = self.lla_to_ecef(ref_llacoordinates)

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
        https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_geodetic_to_ECEF_coordinates
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

        n = a / sqrt(1 - e ** 2 * sin(latitude) ** 2)

        x = (n + altitude) * cos(latitude) * cos(longitude)
        y = (n + altitude) * cos(latitude) * sin(longitude)
        z = (b ** 2 / a ** 2 * n + altitude) * sin(latitude)

        return x, y, z

    def ecef_to_lla(self, coordinates):
        """Convert from ECEF coordinates to LLA coordinates

        ECEF: Earth-Centered, Earth-Fixed
        LLA: Latitude, Longitude, Altitude

        The conversion formulas are taken from
        https://gist.github.com/klucar/1536054

        :param coordinates: tuple of X, Y, and Z (in meters).
        :return: latitude, longitude (in degrees) and altitude (in meters).

        """
        x, y, z = coordinates

        a = self.geode.a
        b = self.geode.b
        e = self.geode.e
        eprime = self.geode.eprime

        p = sqrt(x ** 2 + y ** 2)
        th = atan2(a * z, b * p)

        longitude = atan2(y, x)
        latitude = atan2((z + eprime ** 2 * b * sin(th) ** 3),
                         (p - e ** 2 * a * cos(th) ** 3))
        n = a / sqrt(1 - e ** 2 * sin(latitude) ** 2)
        altitude = p / cos(latitude) - n

        return degrees(latitude), degrees(longitude), altitude

    def ecef_to_enu(self, coordinates):
        """Convert from ECEF coordinates to ENU coordinates

        ECEF: Earth-Centered, Earth-Fixed
        ENU: East, North, Up

        The conversion formulas are taken from
        https://en.wikipedia.org/wiki/Geodetic_system#From_ECEF_to_ENU

        :param coordinates: a tuple containing the ECEF coordinates (in meters)
                            of the point to transform
        :return: east, north, and up (in meters).

        """
        latitude, longitude, altitude = self.ref_lla
        xr, yr, zr = self.ref_ecef
        x, y, z = coordinates

        lat = radians(latitude)
        lon = radians(longitude)

        transformation = array([
            [           -sin(lon),             cos(lon),       0.],  # noqa
            [-sin(lat) * cos(lon), -sin(lat) * sin(lon), cos(lat)],
            [ cos(lat) * cos(lon),  cos(lat) * sin(lon), sin(lat)]])  # noqa

        coordinates = array([x - xr, y - yr, z - zr])

        return tuple(transformation.dot(coordinates))

    def enu_to_ecef(self, coordinates):
        """Convert from ENU coordinates to ECEF coordinates

        ENU: East, North, Up
        ECEF: Earth-Centered, Earth-Fixed

        :param coordinates: a tuple containing the ENU coordinates (in meters).
        :return: ECEF coordinates (in meters).

        """
        latitude, longitude, altitude = self.ref_lla
        xr, yr, zr = self.ref_ecef

        lat = radians(latitude)
        lon = radians(longitude)

        transformation = array([
            [-sin(lon), -sin(lat) * cos(lon), cos(lat) * cos(lon)],
            [ cos(lon), -sin(lat) * sin(lon), cos(lat) * sin(lon)],  # noqa
            [       0.,             cos(lat),            sin(lat)]])  # noqa

        x, y, z = transformation.dot(array(coordinates))

        return x + xr, y + yr, z + zr

    def __repr__(self):
        return f"{self.__class__.__name__}({self.ref_lla!r})"
