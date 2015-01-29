""" HiSPARC cluster definitions

    The BaseCluster class defines a HiSPARC cluster consisting of one or
    more stations.  The Station class defines a HiSPARC station,
    consisting of one or more detectors.

"""
from __future__ import division

from math import sqrt, pi, sin, cos, atan2, radians
import warnings

from numpy import mean

from .transformations import axes, geographic
from . import api


class Detector(object):
    """A HiSPARC detector"""

    _detector_size = (.5, 1.)

    def __init__(self, station, position, orientation='UD'):
        """Initialize detector

        :param station: station instance this detector is part of
        :param position: x,y,z position of the center of the detectors
                         relative to the station center. z is optional.
        :param orientation: either 'UD' or 'LR' meaning an up-down or
                            left-right orientation of the long side of
                            the detector respectively.

        """
        self.station = station
        self.x = position[0]
        self.y = position[1]
        self.z = position[2] if len(position) == 3 else 0.
        self.orientation = orientation

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

        x = X + (self.x * cosa - self.y * sina)
        y = Y + (self.x * sina + self.y * cosa)
        z = Z + self.z

        return x, y, z

    def get_polar_coordinates(self):
        r, phi, _ = self.get_cylindrical_coordinates()
        return r, phi

    def get_cylindrical_coordinates(self):
        x, y, z = self.get_coordinates()
        r, phi, z = axes.cartesian_to_cylindrical(x, y, z)
        return r, phi, z

    def get_corners(self):
        """Get the x, y coordinates of the detector corners

        The z-coordinate is not returned because all detectors are
        assumed to be laying flat.

        :return: x, y coordinates of detector corners
        :rtype: list of (x, y) tuples

        """
        X, Y, _, alpha = self.station.get_coordinates()

        x = self.x
        y = self.y
        orientation = self.orientation
        size = self.detector_size

        dx = size[0] / 2
        dy = size[1] / 2

        if orientation == 'UD':
            corners = [(x - dx, y - dy), (x + dx, y - dy),
                       (x + dx, y + dy), (x - dx, y + dy)]
        elif orientation == 'LR':
            corners = [(x - dy, y - dx), (x + dy, y - dx),
                       (x + dy, y + dx), (x - dy, y + dx)]
        else:
            raise Exception("Unknown detector orientation: %s" % orientation)

        sina = sin(alpha)
        cosa = cos(alpha)
        corners = [[x * cosa - y * sina, x * sina + y * cosa] for x, y in
                   corners]

        return [(X + x, Y + y) for x, y in corners]


class Station(object):
    """A HiSPARC station"""

    _detectors = None

    def __init__(self, cluster, station_id, position, angle=0, detectors=None,
                 number=None):
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
        self.x = position[0]
        self.y = position[1]
        self.z = position[2] if len(position) == 3 else 0.
        self.angle = angle
        self.number = number if number else station_id

        if detectors is None:
            # detector positions for a standard station
            station_size = 10
            a = station_size / 2
            b = a * sqrt(3)
            detectors = [((0, b, 0), 'UD'), ((0, b / 3, 0), 'UD'),
                         ((-a, 0, 0), 'LR'), ((a, 0, 0), 'LR')]

        for position, orientation in detectors:
            self._add_detector(position, orientation)

    def _add_detector(self, position, orientation):
        """Add detector to station

        :param position: x,y,z positions of the center of the detectors
                         relative to the station center.
        :param orientation: Orientation is either 'UD' or 'LR' meaning
                            an up-down or left-right orientation of the
                            long side of the detector respectively.

        """
        if self._detectors is None:
            self._detectors = []
        self._detectors.append(Detector(self, position, orientation))

    @property
    def detectors(self):
        return self._detectors

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

        x = X + (self.x * cosa - self.y * sina)
        y = Y + (self.x * sina + self.y * cosa)
        z = Z + self.z
        alpha = alpha + self.angle

        return x, y, z, alpha

    def get_polar_alpha_coordinates(self):
        r, phi, _, alpha = self.get_cylindrical_alpha_coordinates()
        return r, phi, alpha

    def get_cylindrical_alpha_coordinates(self):
        x, y, z, alpha = self.get_coordinates()
        r, phi, z = axes.cartesian_to_cylindrical(x, y, z)
        return r, phi, z, alpha

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

        x0 = mean(x)
        y0 = mean(y)
        z0 = mean(z)

        return x0, y0, z0


class BaseCluster(object):
    """Base class for HiSPARC clusters"""

    _stations = None

    def __init__(self, position=(0, 0, 0), angle=0):
        """Override this function to build your cluster

        :param position: x,y,z position for the center of the cluster
        :param angle: rotation of the cluster in the x,y-plane.

        """
        self.x = position[0]
        self.y = position[1]
        self.z = position[2] if len(position) == 3 else 0.
        self.alpha = angle

    def _add_station(self, position, angle, detectors=None, number=None):
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
                                      detectors, number))

    @property
    def stations(self):
        return self._stations

    def get_station(self, number):
        """Get a station by its number"""

        for station in self._stations:
            if number == station.number:
                return station

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

        :param s0, s1: The station ids for the two stations.
        :return: r, phi, z; the distance between and direction of the two
            given stations.

        """
        x0, y0, z0, alpha0 = self.stations[s0].get_coordinates()
        x1, y1, z1, alpha1 = self.stations[s1].get_coordinates()

        r = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        phi = atan2((y1 - y0), (x1 - x0))
        z = z1 - z0

        return r, phi, z

    def calc_xy_center_of_mass_coordinates(self):
        x, y, _ = self.calc_center_of_mass_coordinates()
        return x, y

    def calc_center_of_mass_coordinates(self):
        """Calculate center of mass coordinates of all detectors in cluster

        :return: x, y; coordinates of cluster center relative to
            absolute coordinate system

        """
        x, y, z = zip(*[detector.get_coordinates()
                      for station in self.stations
                      for detector in station.detectors])

        x0 = mean(x)
        y0 = mean(y)
        z0 = mean(z)

        return x0, y0, z0


class RAlphaBetaStations(BaseCluster):

    """Add detectors to stations using r, alpha coordinates

    This is meant for data from the Publicdb Database API, which
    uses that coordinate system.

    """

    def _add_station(self, position, detectors, number=None):
        """Add a station to the cluster

        :param x,y,z: position of the station relative to cluster center
        :param detectors: list of tuples. Each tuple consists of (r, alpha,
            z, beta) these define the position of the detector relative to
            the GPS and the orientation of the detector. r is the
            distance between the center of the scintillator and the GPS
            in meters in the x,y-plane. alpha is the clock-wise turning
            angle between North and the detector relative to the GPS in
            degrees. z is the height of the detector relative to the
            GPS. beta is, like alpha, the clock-wise turning angle
            between North and the long side of the detector.
        :param number: optional unique identifier for a station this can
            later be used to find a specific station and makes it easier
            to link to a real station. If not given it will be equal to
            the station_id generated by this function. Either use this
            for all added stations or for none.

        Beta for detectors is currently ignored.

        Example::

            >>> cluster = RAlphaBetaStations()
            >>> cluster._add_station((0, 0, 0), [(5, 90, 0, 0), (5, 270, 0, 0)], 104)

        """
        detectors = [((sin(radians(alpha)) * r, cos(radians(alpha)) * r, z),
                      'UD') for r, alpha, z, beta in detectors]

        super(RAlphaBetaStations, self)._add_station(position, 0, detectors,
                                                     number)


class SimpleCluster(BaseCluster):
    """Define a simple cluster containing four stations"""

    def __init__(self, size=250):
        """Build the cluster"""

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
        """Build the cluster"""

        super(SingleStation, self).__init__()

        self._add_station((0, 0, 0), 0)


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


class ScienceParkCluster(BaseCluster):

    """A cluster containing stations from the Science Park subcluster

    :param stations: A list of station numbers to include. Only stations
        from the Science Park subcluster are supported. By default 507
        and 509 are excluded.

    """

    def __init__(self, stations=[501, 502, 503, 504, 505, 506, 508]):
        super(ScienceParkCluster, self).__init__()

        try:
            gps_coordinates = {}
            for station in stations:
                coordinates = api.Station(station).location()
                gps_coordinates[station] = (coordinates['latitude'],
                                            coordinates['longitude'],
                                            coordinates['altitude'])
        except:
            warnings.warn('Could not get values from the server, Using '
                          'hard-coded values.', UserWarning)
            # 1 day self-survey (8 april 2011) + 506 (Niels, pos from site on
            # 2 dec, 2011) + 508/509 (from site on 8 jul 2013)
            gps_coordinates = {
                501: (52.355924173294305, 4.951144021644267, 56.1023459415882),
                502: (52.355293344895919, 4.9501047083812697, 55.954367009922),
                503: (52.356254735127557, 4.9529437445598328, 51.582641703076),
                504: (52.357178777910278, 4.9543838852175561, 54.622688433155),
                505: (52.357251580629246, 4.9484007564706891, 47.730995402671),
                506: (52.3571787512, 4.95198605591, 43.8700314863),
                507: (52.3560055099, 4.95147879159, 56.7735242238),
                508: (52.3563513341, 4.95070840124, 52.51091104),
                509: (52.3545582682, 4.95569730394, 59.942809986)}

        # 502, 505, 508 are now diamond shapes, rotation has less
        # meaning, need positions of every detector to GPS
        station_rotations = {501: 135, 502: -15, 503: 45, 504: 175, 505: 86,
                             506: 267, 507: 0, 508: -135, 509: 135}

        reference = gps_coordinates[stations[0]]
        transformation = geographic.FromWGS84ToENUTransformation(reference)

        for station in stations:
            enu = transformation.transform(gps_coordinates[station])
            alpha = station_rotations[station] / 180 * pi

            if station not in [501, 502, 505, 508]:
                detectors = [((0., 8.66), 'UD'), ((0., 2.89), 'UD'),
                             ((-5., 0.), 'LR'), ((5., 0.), 'LR')]
            elif station == 501:
                # Precise position measurement of 501
                detectors = [((0.37, 8.62), 'UD'), ((.07, 2.15), 'UD'),
                             ((-5.23, 0.), 'LR'), ((5.08, 0.), 'LR')]
            elif station == 502:
                # 502 is (since 17 October 2011) diamond-shaped,
                # with detector 2 moved to the side in LR orientation.
                # Furthermore, detectors 3 and 4 are reversed (cabling issue)
                station_size = 10.
                a = station_size / 2
                b = a * sqrt(3)
                detectors = [((0., b), 'UD'), ((a * 2, b), 'LR'),
                             ((a, 0.), 'LR'), ((-a, 0.), 'LR')]
            elif station == 505:
                # 505 is (since 24 April 2013) square-shaped,
                # detector 1 is moved to the left and detector 2 next to it.
                station_size = 10.
                a = station_size / 2
                detectors = [((-a, station_size), 'UD'), ((a, station_size), 'UD'),
                             ((-a, 0.), 'LR'), ((a, 0.), 'LR')]
            elif station == 508:
                # 508 is diamond-shaped, with detector 2 moved to the
                # side of detector 1 in UD orientation.
                station_size = 10.
                a = station_size / 2
                b = a * sqrt(3)
                detectors = [((0., b), 'UD'), ((a * 2, b), 'UD'),
                             ((-a, 0.), 'LR'), ((a, 0.), 'LR')]
            else:
                raise RuntimeError("Programming error. Station unknown.")

            self._add_station(enu, alpha, detectors, station)


class HiSPARCStations(RAlphaBetaStations):

    """A cluster containing any real station from the HiSPARC network

    The gps position and number of detectors are taken from the API.
    The detector positions are retrieved if available, otherwise
    default values are used!

    :param stations: A list of station numbers to include. The
        coordinates are retrieved from the Public Database API.
        The first station is placed at the origin of the cluster.

    Example::

        >>> cluster = HiSPARCStations([7001, 7002, 7003])

    """

    def __init__(self, stations, allow_missing=False):
        super(HiSPARCStations, self).__init__()

        for i, station in enumerate(stations):
            try:
                station_info = api.Station(station)
            except:
                if allow_missing:
                    warnings.warn('Could not get info for station %d, '
                                  'using GPS location (0, 0, 0).' % station,
                                  UserWarning)
                    n_detectors = 2
                    detectors = n_detectors * [{'radius': 0, 'alpha': 0,
                                                'height': 0, 'beta': 0}]
                    gps = (0., 0., 0.)
                else:
                    raise KeyError('Could not get info for station %d.' %
                                   station)
            else:
                location = station_info.location()
                n_detectors = station_info.n_detectors()
                detectors = station_info.detectors()
                gps = (location['latitude'],
                       location['longitude'],
                       location['altitude'])

            if i == 0:
                transformation = geographic.FromWGS84ToENUTransformation(gps)

            enu = transformation.transform(gps)

            # Default positions in (r, alpha, z, beta)
            if n_detectors == 2:
                default_detectors = [(5, 90, 0, 0), (5, 270, 0, 0)]
            elif n_detectors == 4:
                default_detectors = [(8.66, 0, 0, 0), (2.89, 0, 0, 0),
                                     (5, -90, 0, 90), (5, 90, 0, 90)]
            else:
                raise RuntimeError("Detector count unknown for station %d." %
                                   station)
            detectors = [
                (d['radius'] if d['radius'] else default_detectors[i][0],
                 d['alpha'] if d['alpha'] else default_detectors[i][1],
                 d['height'] if d['height'] else default_detectors[i][2],
                 d['beta'] if d['beta'] else default_detectors[i][3])
                for i, d in enumerate(detectors)]

            self._add_station(enu, detectors, station)

        warnings.warn('Detector positions may be wrong, defaults are used '
                      'when not available from API!', UserWarning)


class HiSPARCNetwork(HiSPARCStations):

    """A cluster containing all station from the HiSPARC network"""

    def __init__(self):
        network = api.Network()
        stations = [station for station in network.station_numbers()]
        super(HiSPARCNetwork, self).__init__(stations, allow_missing=True)
