""" HiSPARC cluster definitions

    The BaseCluster class defines a HiSPARC cluster consisting of one or
    more stations.  The Station class defines a HiSPARC station,
    consisting of one or more detectors.

"""
from __future__ import division

from math import sqrt, pi


class Station(object):
    """A HiSPARC station"""

    __detector_size = (.5, 1.)

    def __init__(self, position, angle, detectors):
        self.position = position
        self.angle = angle
        self.detectors = detectors

    @property
    def detector_size(self):
        return self.__detector_size


class BaseCluster(object):
    """Base class for HiSPARC clusters"""

    __stations = None

    def __init__(self):
        """Override this function to build your cluster"""
        pass

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
        self.__stations.append(Station(position, angle, detectors))

    @property
    def stations(self):
        return self.__stations


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
        self._add_station(( 0,  0),      0,          detectors)
        self._add_station(( 0,  2 * B),  0,          detectors)
        self._add_station((-A, -B),      2 * pi / 3, detectors)
        self._add_station(( A, -B),     -2 * pi / 3, detectors)
