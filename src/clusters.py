""" HiSPARC cluster definitions

    The BaseCluster class defines a HiSPARC cluster consisting of one or
    more stations.  The Station class defines a HiSPARC station,
    consisting of one or more detectors.

"""
from __future__ import division

from math import sqrt, pi, sin, cos


class Detector(object):
    """A HiSPARC detector"""

    __detector_size = (.5, 1.)

    def __init__(self, station, x, y, orientation):
        """Initialize detector

        :param station: station instance this detector is part of
        :param x, y, orientation: x and y are x and y positions of the
            center of the detectors relative to the station center.
            Orientation is either 'UD' or 'LR' meaning an up-down or
            left-right orientation of the long side of the detector
            respectively.

        """
        self.station = station
        self.x = x
        self.y = y
        self.orientation = orientation

    @property
    def detector_size(self):
        return self.__detector_size

    def get_corners(self, X, Y, alpha):
        """Get the x, y coordinates of the detector corners

        :param X, Y: X, Y coordinates of center of station
        :param alpha: rotation angle of entire station

        :return: x, y coordinates of detector corners
        :rtype: list of (x, y) tuples

        """
        x = self.x
        y = self.y
        orientation = self.orientation
        size = self.detector_size

        dx = size[0] / 2
        dy = size[1] / 2

        if orientation == 'UD':
            corners = [(x - dx, y - dy), (x + dx, y - dy), (x + dx, y + dy),
                       (x - dx, y + dy)]
        elif orientation == 'LR':
            corners = [(x - dy, y - dx), (x + dy, y - dx), (x + dy, y + dx),
                       (x - dy, y + dx)]
        else:
            raise Exception("Unknown detector orientation: %s" % orientation)

        if alpha is not None:
            sina = sin(alpha)
            cosa = cos(alpha)
            corners = [[x * cosa - y * sina, x * sina + y * cosa] for x, y in
                       corners]

        return [(X + x, Y + y) for x, y in corners]


class Station(object):
    """A HiSPARC station"""

    __detectors = None

    def __init__(self, cluster, position, angle, detectors):
        """Initialize station

        :param cluster: cluster this station is a part of
        :param position: tuple of (x, y) values
        :param angle: angle of rotation of the station in radians
        :param detectors: list of tuples.  Each tuple consists of (dx, dy,
            orientation) where dx and dy are x and y positions of the
            center of the detectors relative to the station center.
            Orientation is either 'UD' or 'LR' meaning an up-down or
            left-right orientation of the long side of the detector
            respectively.

        """
        self.cluster = cluster
        self.position = position
        self.angle = angle
        for x, y, orientation in detectors:
            self._add_detector(x, y, orientation)

    def _add_detector(self, x, y, orientation):
        """Add detector to station

        :param x, y, orientation: x and y are x and y positions of the
            center of the detectors relative to the station center.
            Orientation is either 'UD' or 'LR' meaning an up-down or
            left-right orientation of the long side of the detector
            respectively.

        """
        if self.__detectors is None:
            self.__detectors = []
        self.__detectors.append(Detector(self, x, y, orientation))

    @property
    def detectors(self):
        return self.__detectors

    def get_xyalpha_coordinates(self):
        """Calculate coordinates of a station

        :return: x, y, alpha; coordinates and rotation of station relative to
            absolute coordinate system

        """
        X, Y, alpha = self.cluster.get_xyalpha_coordinates()

        sx, sy = self.position
        xp = sx * cos(alpha) - sy * sin(alpha)
        yp = sx * sin(alpha) + sy * cos(alpha)

        x = X + xp
        y = Y + yp
        angle = alpha + self.angle

        return x, y, angle


class BaseCluster(object):
    """Base class for HiSPARC clusters"""

    __stations = None

    def __init__(self, position=(0., 0.), angle=0.):
        """Override this function to build your cluster"""
        self._x, self._y = position
        self._alpha = angle

    def _add_station(self, position, angle, detectors):
        """Add a station to the cluster

        :param position: tuple of (x, y) values
        :param angle: angle of rotation of the station in radians
        :param detectors: list of tuples.  Each tuple consists of (dx, dy,
            orientation) where dx and dy are x and y positions of the
            center of the detectors relative to the station center.
            Orientation is either 'UD' or 'LR' meaning an up-down or
            left-right orientation of the long side of the detector
            respectively.

        Example::

            >>> cluster = BaseCluster()
            >>> cluster._add_station((0, 0), pi / 2, [(-5, 0, 'UD'), (5, 0, 'UD')])
        """
        # Need to make __stations an instance variable to be able to
        # pickle it.  An assignment takes care of that.
        if self.__stations is None:
            self.__stations = []
        self.__stations.append(Station(self, position, angle, detectors))

    @property
    def stations(self):
        return self.__stations

    def get_xyalpha_coordinates(self):
        return self._x, self._y, self._alpha

    def set_xyalpha_coordinates(self, position):
        self._x, self._y, self._alpha = position

    def set_rphialpha_coordinates(self, position):
        r, phi, self._alpha = position
        self._x = r * cos(phi)
        self._y = r * sin(phi)


class SimpleCluster(BaseCluster):
    """Define a simple cluster containing four stations"""

    def __init__(self, size=250):
        """Build the cluster"""

        # calculate detector positions for a four-detector station
        station_size = 10
        a = station_size / 2
        b = a / 3 * sqrt(3)
        detectors = [(0., 2 * b, 'UD'), (0., 0., 'UD'),
                     (-a, -b, 'LR'), (a, -b, 'LR')]

        # calculate station positions.  the cluster resembles a single
        # four-detector HiSPARC station, but scaled up
        A = size / 2
        B = A / 3 * sqrt(3)
        self._add_station((0, 2 * B), 0, detectors)
        self._add_station((0, 0), 0, detectors)
        self._add_station((-A, -B), 2 * pi / 3, detectors)
        self._add_station((A, -B), -2 * pi / 3, detectors)
