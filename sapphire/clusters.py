""" Define HiSPARC detectors, stations and clusters.

    The BaseCluster class defines a HiSPARC cluster consisting of one or
    more stations.  The Station class defines a HiSPARC station,
    consisting of one or more Detectors.

"""
from __future__ import division

from math import sqrt, pi, sin, cos, atan2
import warnings

import numpy as np

from .transformations import axes, geographic
from . import api
from .utils import get_active_index


class Detector(object):
    """A HiSPARC detector"""

    _detector_size = (.5, 1.)

    def __init__(self, station, position, orientation='UD',
                 detector_timestamps=[0]):
        """Initialize detector

        :param station: station instance this detector is part of.
        :param position: x,y,z position of the center of the detectors
            relative to the station center. z is optional.
        :param orientation: orientation of the long side of the detector.
            Either the angle in radians, or 'UD' or 'LR' meaning an
            up-down or left-right orientation of the long side of the
            detector respectively.

        """
        self.station = station
        if hasattr(position[0], "__len__"):
            self.x = position[0]
            self.y = position[1]
            self.z = position[2] if len(position) == 3 else [0.] * len(self.x)
        else:
            self.x = [position[0]]
            self.y = [position[1]]
            self.z = [position[2]] if len(position) == 3 else [0.]
        if orientation is 'UD':
            self.orientation = [0] * len(self.x)
        elif orientation is 'LR':
            self.orientation = [pi / 2] * len(self.x)
        else:
            if hasattr(orientation, "__len__"):
                self.orientation = orientation
            else:
                self.orientation = [orientation]
        if len(detector_timestamps) == len(self.x):
            self.timestamps = detector_timestamps
        else:
            raise Exception('Number of timestamps must equal number of '
                            'postions')
        self.index = -1

    def _update_timestamp(self, timestamp):
        """Get the position index valid for the given timestamp

        :param timestamp: timestamp in seconds.

        """
        self.index = get_active_index(self.timestamps, timestamp)

    @property
    def detector_size(self):
        return self._detector_size

    def get_area(self):
        return self._detector_size[0] * self._detector_size[1]

    def get_xy_coordinates(self):
        x, y, _ = self.get_coordinates()
        return x, y

    def get_coordinates(self):
        X, Y, Z, alpha = self.station.get_coordinates()

        sina = sin(alpha)
        cosa = cos(alpha)

        x = X + (self.x[self.index] * cosa - self.y[self.index] * sina)
        y = Y + (self.x[self.index] * sina + self.y[self.index] * cosa)
        z = Z + self.z[self.index]

        return x, y, z

    def get_polar_coordinates(self):
        r, phi, _ = self.get_cylindrical_coordinates()
        return r, phi

    def get_cylindrical_coordinates(self):
        x, y, z = self.get_coordinates()
        r, phi, z = axes.cartesian_to_cylindrical(x, y, z)
        return r, phi, z

    def get_lla_coordinates(self):
        """Get detector LLA coordinates (latitude, longitude, altitude)

        get detector WGS84 coordinates
        LLA: Latitude, Longitude, Altitude

        :return: tuple (latitude, longitude, altitude).
                 Latitude, longitude in degrees. Altitude in meters.

        """
        lla = self.station.cluster.lla
        enu = self.get_coordinates()

        transform = geographic.FromWGS84ToENUTransformation(lla)
        latitude, longitude, altitude = transform.enu_to_lla(enu)

        return latitude, longitude, altitude

    def get_corners(self):
        """Get the x, y coordinates of the detector corners

        The z-coordinate is not returned because all detectors are
        assumed to be laying flat.

        :return: coordinates of detector corners, list of (x, y) tuples.

        """
        X, Y, _, alpha = self.station.get_coordinates()

        x = self.x[self.index]
        y = self.y[self.index]
        o = self.orientation[self.index]
        size = self.detector_size

        # detector frame
        dx = size[0] / 2
        dy = size[1] / 2
        corners = [(-dx, -dy), (dx, -dy), (dx, dy), (-dx, dy)]

        # station frame
        coso = cos(-o)
        sino = sin(-o)
        corners = [(x + cx * coso - cy * sino, y + cx * sino + cy * coso)
                   for cx, cy in corners]

        # cluster frame
        sina = sin(alpha)
        cosa = cos(alpha)
        corners = [(X + xc * cosa - yc * sina, Y + xc * sina + yc * cosa)
                   for xc, yc in corners]

        return corners


class Station(object):
    """A HiSPARC station"""

    _detectors = None

    def __init__(self, cluster, station_id, position, angle=None,
                 detectors=None, station_timestamps=[0],
                 detector_timestamps=[0], number=None):
        """Initialize station

        :param cluster: cluster this station is a part of
        :param station_id: int (unique identifier)
        :param position: x,y,z position of the station center relative
            to the cluster center, z is optional.
        :param angle: angle of rotation of the station in radians
        :param detectors: list of tuples.  Each tuple consists of (dx, dy,
            orientation) where dx and dy are x and y positions of the
            center of the detectors relative to the station center.
            Orientation is either 'UD' or 'LR' meaning an up-down or
            left-right orientation of the long side of the detector
            respectively.
        :param number: optional unique identifier for a station this can
            be used by the cluster to find a specific station and makes
            it easier to link to a real station. If not given it will be
            equal to the station_id.

        """
        self.cluster = cluster
        self.station_id = station_id
        if hasattr(position[0], "__len__"):
            self.x = position[0]
            self.y = position[1]
            self.z = position[2] if len(position) == 3 else [0.] * len(self.x)
        else:
            self.x = [position[0]]
            self.y = [position[1]]
            self.z = [position[2]] if len(position) == 3 else [0.]
        if angle is None:
            self.angle = [0.] * len(self.x)
        elif hasattr(angle, "__len__"):
            self.angle = angle
        else:
            self.angle = [angle]
        self.number = number if number else station_id

        if len(station_timestamps) == len(self.x):
            self.timestamps = station_timestamps
        else:
            raise Exception('Number of timestamps must equal number of '
                            'postions')

        if detectors is None:
            # detector positions for a standard station
            station_size = 10
            a = station_size / 2
            b = a * sqrt(3)
            detectors = [((0, b, 0), 'UD'), ((0, b / 3, 0), 'UD'),
                         ((-a, 0, 0), 'LR'), ((a, 0, 0), 'LR')]

        for position, orientation in detectors:
            self._add_detector(position, orientation, detector_timestamps)
        self.index = -1

    def _update_timestamp(self, timestamp):
        """Get the position index valid for the given timestamp

        :param timestamp: timestamp in seconds.

        """
        for detector in self.detectors:
            detector._update_timestamp(timestamp)
        self.index = get_active_index(self.timestamps, timestamp)

    def _add_detector(self, position, orientation, detector_timestamps):
        """Add detector to station

        :param position: x,y,z positions of the center of the detectors
                         relative to the station center.
        :param orientation: Orientation is either 'UD' or 'LR' meaning
                            an up-down or left-right orientation of the
                            long side of the detector respectively.

        """
        if self._detectors is None:
            self._detectors = []
        self._detectors.append(Detector(self, position, orientation,
                                        detector_timestamps))

    @property
    def detectors(self):
        return self._detectors

    def get_area(self, detector_ids=None):
        """Get the total area covered by the detectors

        :param detector_ids: list of detectors for which to get the total area.
        :return: total area of the detectors in m^2.

        """
        if detector_ids is not None:
            return sum(self._detectors[id].get_area() for id in detector_ids)
        else:
            return sum(d.get_area() for d in self._detectors)

    def get_xy_coordinates(self):
        """Same as get_coordinates but without the z and alpha"""
        x, y, _, _ = self.get_coordinates()
        return x, y

    def get_xyalpha_coordinates(self):
        """Same as get_coordinates but without the z"""
        x, y, _, alpha = self.get_coordinates()
        return x, y, alpha

    def get_coordinates(self):
        """Calculate coordinates of a station

        :return: x, y, z, alpha; coordinates and rotation of station
                 relative to absolute coordinate system

        """
        X, Y, Z, alpha = self.cluster.get_coordinates()

        sina = sin(alpha)
        cosa = cos(alpha)

        x = X + (self.x[self.index] * cosa - self.y[self.index] * sina)
        y = Y + (self.x[self.index] * sina + self.y[self.index] * cosa)
        z = Z + self.z[self.index]
        alpha = alpha + self.angle[self.index]

        return x, y, z, alpha

    def get_polar_alpha_coordinates(self):
        r, phi, _, alpha = self.get_cylindrical_alpha_coordinates()
        return r, phi, alpha

    def get_cylindrical_alpha_coordinates(self):
        x, y, z, alpha = self.get_coordinates()
        r, phi, z = axes.cartesian_to_cylindrical(x, y, z)
        return r, phi, z, alpha

    def get_lla_coordinates(self):
        """Get station LLA coordinates (latitude, longitude, altitude)

        LLA: Latitude, Longitude, Altitude

        :return: tuple (latitude, longitude, altitude).
                 Latitude, longitude in degrees. Altitude in meters.

        """
        lla = self.cluster.lla
        x, y, z, alpha = self.get_coordinates()
        enu = (x, y, z)

        transform = geographic.FromWGS84ToENUTransformation(lla)
        latitude, longitude, altitude = transform.enu_to_lla(enu)
        latitude = latitude if abs(latitude) > 1e-7 else 0.
        longitude = longitude if abs(longitude) > 1e-7 else 0.
        altitude = altitude if abs(altitude) > 1e-7 else 0.

        return latitude, longitude, altitude

    def calc_r_and_phi_for_detectors(self, d0, d1):
        r, phi, _ = self.calc_rphiz_for_detectors(d0, d1)
        return r, phi

    def calc_rphiz_for_detectors(self, d0, d1):
        """Calculate angle and distance between detectors

        :param d0,d1: detector ids to find the vector between.
        :return: r,phi,z pointing from d0 to d1.

        """
        x0, y0, z0 = self.detectors[d0].get_coordinates()
        x1, y1, z1 = self.detectors[d1].get_coordinates()

        r = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        phi = atan2((y1 - y0), (x1 - x0))
        dz = z1 - z0

        return r, phi, dz

    def calc_xy_center_of_mass_coordinates(self):
        x, y, _ = self.calc_center_of_mass_coordinates()
        return x, y

    def calc_center_of_mass_coordinates(self):
        """Calculate center of mass coordinates of detectors in station

        :return: x, y, z; coordinates of station center relative to
            absolute coordinate system

        """
        x, y, z = zip(*[detector.get_coordinates()
                      for detector in self.detectors])

        x0 = np.nanmean(x)
        y0 = np.nanmean(y)
        z0 = np.nanmean(z)

        return x0, y0, z0


class BaseCluster(object):
    """Base class for HiSPARC clusters"""

    _stations = None

    def __init__(self, position=(0, 0, 0), angle=0,
                 lla=(52.35592417, 4.95114402, 56.10234594)):
        """Override this function to build your cluster

        :param position: x,y,z position for the center of the cluster.
        :param angle: rotation of the cluster in the x,y-plane.
        :param lla: Reference WGS84 location of the cluster origin.
                    Defaults to (old) GPS location of station 501 (Nikhef).

        """
        self.x = position[0]
        self.y = position[1]
        self.z = position[2] if len(position) == 3 else 0.
        self.alpha = angle
        self.lla = lla
        # Set initial timestamp in the future to use latest positions
        # 2 ** 31 - 1 == 19 Jan 2038
        self._timestamp = 2147483647

    def set_timestamp(self, timestamp):
        """Set the timestamp to set the active station and detector locations

        :param timestamp: timestamp in seconds.

        """
        self._timestamp = timestamp
        for station in self.stations:
            station._update_timestamp(self._timestamp)

    def _add_station(self, position, angle=None, detectors=None,
                     station_timestamps=[0], detector_timestamps=[0],
                     number=None):
        """Add a station to the cluster

        :param position: x,y,z position of the station relative to
            cluster center. z is optional.
        :param angle: angle of rotation of the station in radians
        :param detectors: list of tuples.  Each tuple consists of (dx, dy,
            dz, orientation) where dx, dy and dz are the positions of the
            center of the detectors relative to the station center.
            Orientation is either 'UD' or 'LR' meaning an up-down or
            left-right orientation of the long side of the detector
            respectively.
        :param number: optional unique identifier for a station this can
            later be used to find a specific station and makes it easier
            to link to a real station. If not given it will be equal to
            the station_id generated by this function. Either use this
            for all added stations or for none.

        Example::

            >>> cluster = BaseCluster()
            >>> cluster._add_station((0, 0, 0), pi / 2,
            ...                      [((-5, 0, 0), 'UD'), ((5, 0, 0), 'UD')])

        """
        # Need to make _stations an instance variable to be able to
        # pickle it.  An assignment takes care of that.
        if self._stations is None:
            self._stations = []

        station_id = len(self._stations)
        self._stations.append(Station(self, station_id, position, angle,
                                      detectors, station_timestamps,
                                      detector_timestamps, number))

    @property
    def stations(self):
        return self._stations

    def get_station(self, number):
        """Get a station by its number"""

        for station in self._stations:
            if number == station.number:
                return station

    def get_xy_coordinates(self):
        """Same as get_coordinates but without the z and alpha"""
        return self.x, self.y

    def get_xyalpha_coordinates(self):
        """Like get_coordinates, but without z"""
        return self.x, self.y, self.alpha

    def get_coordinates(self):
        """Get cluster coordinates (x, y, z, alpha).

        The coordinates should be interpreted as follows: first, the
        cluster is rotated over angle alpha, around its original center.
        Then, the cluster is translated to (x, y, z).

        """
        return self.x, self.y, self.z, self.alpha

    def get_polar_alpha_coordinates(self):
        """Like get_cylindrical_coordinates but without z."""
        r, phi, _, alpha = self.get_cylindrical_alpha_coordinates()
        return r, phi, alpha

    def get_cylindrical_alpha_coordinates(self):
        """Get cluster coordinates (r, phi, z, alpha).

        The coordinates should be interpreted as follows: first, the
        cluster is rotated over angle alpha, around its original center.
        Then, the cluster is translated to (r, phi, z).

        """
        r, phi, z = axes.cartesian_to_cylindrical(self.x, self.y, self.z)
        return r, phi, z, self.alpha

    def get_lla_coordinates(self):
        """Get cluster LLA coordinates (latitude, longitude, altitude)

        get cluster WGS84 coordinates
        LLA: Latitude, Longitude, Altitude

        :return: tuple (latitude, longitude, altitude).
                 Latitude, longitude in degrees. Altitude in meters.

        """
        lla = self.lla
        x, y, z, alpha = self.get_coordinates()
        enu = (x, y, z)

        transform = geographic.FromWGS84ToENUTransformation(lla)
        latitude, longitude, altitude = transform.enu_to_lla(enu)

        return latitude, longitude, altitude

    def set_coordinates(self, x, y, z, alpha):
        """Set cluster coordinates (x, y, z, alpha).

        The coordinates should be interpreted as follows: first, the
        cluster is rotated over angle alpha, around its original center.
        Then, the cluster is translated to (x, y, z).

        """
        self.x, self.y, self.z, self.alpha = x, y, z, alpha

    def set_cylindrical_coordinates(self, r, phi, z, alpha):
        """Set cluster coordinates (r, phi, z, alpha).

        The coordinates should be interpreted as follows: first, the
        cluster is rotated over angle alpha, around its original center.
        Than, the cluster is translated to (r, phi).

        """
        self.x, self.y, self.z = axes.cylindrical_to_cartesian(r, phi, z)
        self.alpha = alpha

    def calc_rphiz_for_stations(self, s0, s1):
        """Calculate distance between and direction of two stations

        The calculated result is between the center of mass coordinates
        of the two stations.

        :param s0,s1: The station ids for the two stations.
        :return: r, phi, z; the distance between and direction of the two
            given stations.

        """
        x0, y0, z0 = self.stations[s0].calc_center_of_mass_coordinates()
        x1, y1, z1 = self.stations[s1].calc_center_of_mass_coordinates()

        r = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        phi = atan2((y1 - y0), (x1 - x0))
        z = z1 - z0

        return r, phi, z

    def calc_xy_center_of_mass_coordinates(self):
        """Like calc_center_of_mass_coordinates, but without z"""

        x, y, _ = self.calc_center_of_mass_coordinates()
        return x, y

    def calc_center_of_mass_coordinates(self):
        """Calculate center of mass coordinates of all detectors in cluster

        :return: x, y, z; coordinates of cluster center relative to
            absolute coordinate system

        """
        x, y, z = zip(*[detector.get_coordinates()
                      for station in self.stations
                      for detector in station.detectors])

        x0 = np.mean(x)
        y0 = np.mean(y)
        z0 = np.mean(z)

        return x0, y0, z0

    @staticmethod
    def _distance(c1, c2):
        return np.sqrt(sum((c1 - c2) ** 2))

    def calc_distance_between_stations(self, s1, s2):
        """Calculate distance between two stations

        :param s1,s2: station numbers.
        :return: distance between stations.

        """
        pair = [so for so in self._stations if so.number in (s1, s2)]

        if len(pair) != 2:
            return None

        xyz = [np.array(s.calc_center_of_mass_coordinates()) for s in pair]

        return self._distance(*xyz)

    def calc_horizontal_distance_between_stations(self, s1, s2):
        """Calculate 2D distance between two HiSPARC stations. Ignores altitude

        The 2D plane is the East-North plane defined by the ENU axes at the
        reference location for this cluster.

        :param s1,s2: station numbers.
        :return: distance between stations.

        """
        pair = [so for so in self._stations if so.number in (s1, s2)]

        if len(pair) != 2:
            return None

        xy = [np.array(s.calc_center_of_mass_coordinates()[:-1]) for s in pair]

        return self._distance(*xy)


class CompassStations(BaseCluster):

    """Add detectors to stations using compass coordinates

    Compass coordinates consist of r, alpha, z, beta. These define
    the location and orientation of detectors. For more information
    see `Coordinate systems and units in HiSPARC`.

    This is meant for data from the Publicdb Database API, which
    uses that coordinate system.

    """

    def _add_station(self, position, detectors, station_timestamps=[0],
                     detector_timestamps=[0], number=None):
        """Add a station to the cluster

        :param position: x,y,z coordinates of the station relative
            to cluster center (ENU), can be list of multiple positions.
        :param detectors: list of r,alpha,z,beta coordinates, these
            define the position of the detector relative to the GPS and
            the orientation of the detector. r is the distance between
            the center of the scintillator and the GPS in meters in the
            x,y-plane. alpha is the clock-wise turning angle between
            North and the detector relative to the GPS in degrees. z is
            the height of the detector relative to the GPS. beta is the
            clock-wise turning angle between North and the long side of
            the detector in degrees. Each coordinate may be an array, to
            define multiple station layouts.
        :param station_timestamps: list of timestamps, the timestamp at
            which each of the positions became active.
        :param detector_timestamps: list of timestamps, the timestamp at
            which each of the layouts became active.
        :param number: optional unique identifier for a station this can
            later be used to find a specific station and makes it easier
            to link to a real station. If not given it will be equal to
            the station_id generated by this function. Either use this
            for all added stations or for none.

        Beta for detectors is currently ignored.

        Example::

            >>> cluster = CompassStations()
            >>> cluster._add_station((0, 0, 0), [(7, 0, 1, 0), (7, 90, 0, 0)],
            ...                      number=104)

        """
        detectors = [(axes.compass_to_cartesian(r, alpha, z), np.radians(beta))
                     for r, alpha, z, beta in detectors]

        super(CompassStations, self)._add_station(
            position, None, detectors, station_timestamps, detector_timestamps,
            number)


class SimpleCluster(BaseCluster):

    """Define a simple cluster containing four stations

    :param size: This value is the distance between the three outer stations.

    """

    def __init__(self, size=250):
        super(SimpleCluster, self).__init__()

        # calculate station positions. the cluster resembles a single
        # four-detector HiSPARC station, but scaled up
        A = size / 2
        B = A / sqrt(3)
        self._add_station((0, 2 * B, 0), 0)
        self._add_station((0, 0, 0), 0)
        self._add_station((-A, -B, 0), 2 * pi / 3)
        self._add_station((A, -B, 0), -2 * pi / 3)


class SingleStation(BaseCluster):
    """Define a cluster containing a single station"""

    def __init__(self):
        super(SingleStation, self).__init__()

        self._add_station((0, 0, 0), 0)


class SingleDetectorStation(BaseCluster):
    """Define a cluster containing a single 1-detector station"""

    def __init__(self):
        super(SingleDetectorStation, self).__init__()

        detectors = [((0, 0, 0), 'UD')]

        self._add_station((0, 0, 0), 0, detectors)


class SingleTwoDetectorStation(BaseCluster):
    """Define a cluster containing a single 2 detector station"""

    def __init__(self):
        super(SingleTwoDetectorStation, self).__init__()

        detectors = [((-5, 0, 0), 'UD'), ((5, 0, 0), 'UD')]

        self._add_station((0, 0, 0), 0, detectors)


class SingleDiamondStation(BaseCluster):

    """Define a cluster containing a single diamond shaped station

    Detectors 1, 3 and 4 are in the usual position for a 4 detector
    layout, detector 2 is moved out of the center and positioned to
    create a second equilateral triangle with detectors 1, 2 and 4.

    """

    def __init__(self):
        super(SingleDiamondStation, self).__init__()

        station_size = 10
        a = station_size / 2
        b = a * sqrt(3)
        detectors = [((0., b, 0), 'UD'), ((a * 2, b, 0), 'UD'),
                     ((-a, 0., 0), 'LR'), ((a, 0., 0), 'LR')]

        self._add_station((0, 0, 0), 0, detectors)


class HiSPARCStations(CompassStations):

    """A cluster containing any real station from the HiSPARC network

    The gps position and number of detectors are taken from the API.
    The detector positions are retrieved if available, otherwise
    default values are used!

    :param stations: A list of station numbers to include. The
        coordinates are retrieved from the Public Database API.
        The first station is placed at the origin of the cluster.
    :param skip_missing: Set to True to skip stations which have missing
        location data, otherwise an exception will be raised. Stations
        with missing location data will be excluded. Does not apply
        to missing detector positions.

    Example::

        >>> cluster = HiSPARCStations([7001, 7002, 7003])

    """

    def __init__(self, stations, skip_missing=False, force_fresh=False,
                 force_stale=False):
        super(HiSPARCStations, self).__init__()

        missing_gps = []
        missing_detectors = []

        for i, station in enumerate(stations):
            try:
                station_info = api.Station(station, force_fresh=force_fresh,
                                           force_stale=force_stale)
                locations = station_info.gps_locations
            except:
                if skip_missing:
                    missing_gps.append(station)
                    continue
                else:
                    raise KeyError('Could not get GPS info for station %d.' %
                                   station)
            else:
                llas = locations[['latitude', 'longitude', 'altitude']]
                station_ts = locations['timestamp']
                n_detectors = station_info.n_detectors()

            if i == 0:
                # Most recent location of first station as reference
                self.lla = llas[-1]
                transformation = geographic.FromWGS84ToENUTransformation(
                    self.lla)

            # Station locations in ENU
            enu = [transformation.transform(lla) for lla in llas]
            enu = [list(coordinate) for coordinate in zip(*enu)]

            try:
                detectors = station_info.station_layouts
                fields = ('radius', 'alpha', 'height', 'beta')
                razbs = [[detectors['%s%d' % (field, i)] for field in fields]
                         for i in range(1, n_detectors + 1)]
                detector_ts = detectors['timestamp']
            except:
                missing_detectors.append(station)
                # Fallback detector positions in (r, alpha, z, beta)
                if n_detectors == 2:
                    razbs = [(5, 90, 0, 0), (5, 270, 0, 0)]
                elif n_detectors == 4:
                    d = 10 / sqrt(3)
                    razbs = [(d, 0, 0, 0), (0, 0, 0, 0),
                             (d, -120, 0, 90), (d, 120, 0, 90)]
                else:
                    raise RuntimeError("Detector count unknown for station %d."
                                       % station)
                detector_ts = [0]

            self._add_station(enu, razbs, station_ts, detector_ts, station)

        if len(missing_gps):
            warnings.warn('Could not get GPS location for stations: %s. '
                          'Using (0, 0, 0) instead.' % str(missing_gps))
        if len(missing_detectors):
            warnings.warn('Could not get detector layout for stations %s, '
                          'defaults will be used!' % str(missing_detectors))


class ScienceParkCluster(HiSPARCStations):

    """A cluster containing stations from the Science Park subcluster

    :param stations: A list of station numbers to include. Only stations
        from the Science Park subcluster (5xx) are supported. By default 507
        is excluded.

    """

    def __init__(self, stations=None, skip_missing=False, force_fresh=False,
                 force_stale=False):
        if stations is None:
            network = api.Network(force_fresh, force_stale)
            stations = [sn for sn in network.station_numbers(subcluster=500)
                        if sn != 507]
        else:
            stations = [sn for sn in stations if 500 < sn < 600]
        super(ScienceParkCluster, self).__init__(stations, skip_missing,
                                                 force_fresh, force_stale)


class HiSPARCNetwork(HiSPARCStations):

    """A cluster containing all station from the HiSPARC network"""

    def __init__(self, force_fresh=False, force_stale=False):
        network = api.Network(force_fresh, force_stale)
        stations = network.station_numbers()
        skip_missing = True  # Likely some station without GPS location
        super(HiSPARCNetwork, self).__init__(stations, skip_missing,
                                             force_fresh, force_stale)
