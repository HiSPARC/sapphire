""" Direction reconstruction

    This module contains two classes that can be used to reconstruct
    HiSPARC events and coincidences. The classes know how to extract the
    relevant information from the station and event or cluster and
    coincidence. Various algorithms which do the reconstruction are also
    defined here. The algorithms require positions and arrival times to
    do the reconstruction.

    Each algorithm has a :meth:`~BaseDirectionAlgorithm.reconstruct_common`
    method which always requires arrival times, x, and y positions and
    optionally z positions and previous reconstruction results. The data
    is then prepared for the algorithm and passed to
    the :meth:`~BaseDirectionAlgorithm.reconstruct` method which returns the
    reconstructed theta and phi coordinates.

"""
import warnings
from itertools import combinations

from six.moves import zip_longest
from six import itervalues
from numpy import (nan, isnan, arcsin, arccos, arctan2, sin, cos, tan,
                   sqrt, where, pi, inf, array, cross, dot, sum, zeros)
from scipy.optimize import minimize
from scipy.sparse.csgraph import shortest_path

from .event_utils import (station_arrival_time, detector_arrival_time,
                          relative_detector_arrival_times)
from ..simulations.showerfront import CorsikaStationFront
from ..utils import (pbar, norm_angle, c, make_relative, vector_length,
                     floor_in_base, memoize)
from ..api import Station


NO_OFFSET = [0., 0., 0., 0.]
NO_STATION_OFFSET = (0., 100.)


class EventDirectionReconstruction(object):

    """Reconstruct direction for station events

    This class is aware of 'events' and 'stations'.  Initialize this class
    with a 'station' and you can reconstruct events using
    :meth:`reconstruct_event`. To use other algorithms overwrite the
    ``direct`` and ``fit`` attributes.

    :param station: :class:`sapphire.clusters.Station` object.

    """

    def __init__(self, station):
        self.direct = DirectAlgorithmCartesian3D
        self.fit = RegressionAlgorithm3D
        self.station = station

    def reconstruct_event(self, event, detector_ids=None, offsets=NO_OFFSET,
                          initial=None):
        """Reconstruct a single event

        :param event: an event (e.g. from an events table), or any
            dictionary-like object containing the keys necessary for
            reconstructing the direction of a shower (e.g. arrival times).
        :param detector_ids: list of the detectors to use for
            reconstruction. The detector ids are 0-based, unlike the
            column names in the esd data.
        :param offsets: time offsets for each detector or a
            :class:`~sapphire.api.Station` object.
        :param initial: dictionary with already fitted shower parameters.
        :return: theta, phi, and detector ids.

        """
        t, x, y, z, ids = ([], [], [], [], [])
        if detector_ids is None:
            detector_ids = range(4)
        self.station.cluster.set_timestamp(event['timestamp'])
        if isinstance(offsets, Station):
            offsets = offsets.detector_timing_offset(event['timestamp'])
        for id in detector_ids:
            t_detector = detector_arrival_time(event, id, offsets)
            if not isnan(t_detector):
                dx, dy, dz = self.station.detectors[id].get_coordinates()
                t.append(t_detector)
                x.append(dx)
                y.append(dy)
                z.append(dz)
                ids.append(id)
        if len(t) == 3:
            theta, phi = self.direct.reconstruct_common(t, x, y, z, initial)
        elif len(t) > 3:
            theta, phi = self.fit.reconstruct_common(t, x, y, z, initial)
        else:
            theta, phi = (nan, nan)
        return theta, phi, ids

    def reconstruct_events(self, events, detector_ids=None, offsets=NO_OFFSET,
                           progress=True, initials=None):
        """Reconstruct events

        :param events: the events table for the station from an ESD data file.
        :param detector_ids: detectors to use for the reconstructions.
        :param offsets: time offsets for each detector or a
            :class:`~sapphire.api.Station` object.
        :param progress: if True show a progress bar while reconstructing.
        :param initials: list of dictionaries with already reconstructed shower
                         parameters.
        :return: list of theta, phi, and detector ids.

        """
        if initials is None:
            initials = []
        events = pbar(events, show=progress)
        events_init = zip_longest(events, initials)
        angles = [self.reconstruct_event(event, detector_ids, offsets, initial)
                  for event, initial in events_init]
        if len(angles):
            theta, phi, ids = zip(*angles)
        else:
            theta, phi, ids = ((), (), ())
        return theta, phi, ids

    def __repr__(self):
        return ("<%s, station: %r, direct: %r, fit: %r>" %
                (self.__class__.__name__, self.station, self.direct, self.fit))


class CoincidenceDirectionReconstruction(object):

    """Reconstruct direction for coincidences

    This class is aware of 'coincidences' and 'clusters'.  Initialize
    this class with a 'cluster' and you can reconstruct a coincidence
    using :meth:`reconstruct_coincidence`. To use other algorithms
    overwrite the ``direct``,``fit``, and ``curved`` attributes.

    :param cluster: :class:`~sapphire.clusters.BaseCluster` object.

    """

    def __init__(self, cluster):
        self.direct = DirectAlgorithmCartesian3D
        self.fit = RegressionAlgorithm3D
        self.curved = CurvedRegressionAlgorithm3D()
        self.cluster = cluster

    def reconstruct_coincidence(self, coincidence_events, station_numbers=None,
                                offsets=None, initial=None):
        """Reconstruct a single coincidence

        :param coincidence_events: a coincidence list consisting of three
            or more (station_number, event) tuples.
        :param station_numbers: list of station numbers, to only use
            events from those stations.
        :param offsets: a dictionary of either lists of detector timing
            offsets or :class:`~sapphire.api.Station` objects for each station.
        :param initial: dictionary with already fitted shower parameters.
        :return: list of theta, phi, and station numbers.

        """
        if len(coincidence_events) < 3:
            return nan, nan, []
        if offsets is None:
            offsets = {}
        if initial is None:
            initial = {}

        # Subtract base timestamp to prevent loss of precision
        ts0 = int(coincidence_events[0][1]['timestamp'])
        ets0 = ts0 * int(1e9)
        self.cluster.set_timestamp(ts0)
        t, x, y, z, nums = ([], [], [], [], [])

        offsets = self.get_station_offsets(coincidence_events, station_numbers,
                                           offsets, ts0)

        for station_number, event in coincidence_events:
            if station_numbers is not None:
                if station_number not in station_numbers:
                    continue
            t_off = offsets.get(station_number, NO_OFFSET)
            station = self.cluster.get_station(station_number)
            t_first = station_arrival_time(event, ets0, offsets=t_off,
                                           station=station)
            if not isnan(t_first):
                sx, sy, sz = station.calc_center_of_mass_coordinates()
                t.append(t_first)
                x.append(sx)
                y.append(sy)
                z.append(sz)
                nums.append(station_number)

        if len(t) >= 3 and 'core_x' in initial and 'core_y' in initial:
            theta, phi = self.curved.reconstruct_common(t, x, y, z, initial)
        elif len(t) == 3:
            theta, phi = self.direct.reconstruct_common(t, x, y, z, initial)
        elif len(t) > 3:
            theta, phi = self.fit.reconstruct_common(t, x, y, z, initial)
        else:
            theta, phi = (nan, nan)

        return theta, phi, nums

    def reconstruct_coincidences(self, coincidences, station_numbers=None,
                                 offsets=None, progress=True, initials=None):
        """Reconstruct all coincidences

        :param coincidences: a list of coincidence events, each consisting
                             of three or more (station_number, event) tuples.
        :param station_numbers: list of station numbers, to only use
                                events from those stations.
        :param offsets: dictionary with detector offsets for each station.
                        These detector offsets should be relative to one
                        detector from a specific station.
        :param progress: if True show a progress bar while reconstructing.
        :param initials: list of dictionaries with already reconstructed shower
                         parameters.
        :return: list of theta, phi, and station numbers.

        """
        if offsets is None:
            offsets = {}
        if initials is None:
            initials = []
        coincidences = pbar(coincidences, show=progress)
        coin_init = zip_longest(coincidences, initials)
        angles = [self.reconstruct_coincidence(coincidence, station_numbers,
                                               offsets, initial)
                  for coincidence, initial in coin_init]
        if len(angles):
            theta, phi, nums = zip(*angles)
        else:
            theta, phi, nums = ((), (), ())
        return theta, phi, nums

    def get_station_offsets(self, coincidence_events, station_numbers,
                            offsets, ts0):
        if offsets and isinstance(next(itervalues(offsets)), Station):
            if station_numbers is None:
                # stations in the coincidence
                stations = list({sn for sn, _ in coincidence_events})
            else:
                stations = station_numbers
            midnight_ts = floor_in_base(ts0, 86400)
            offsets = self.determine_best_offsets(stations, midnight_ts,
                                                  offsets)
        return offsets

    @memoize
    def determine_best_offsets(self, station_numbers, midnight_ts, offsets):
        """Determine best combined station and detector offsets

        Check which station is best used as reference. Allow offsets via
        other stations, intermediate stations are used if it reduces the
        offset error.

        :param station_numbers: list of stations in the coincidence or also
                                other stations can are allow to be the
                                reference station.
        :param midnight_ts: timestamp of midnight before the coincidence.
        :param offsets: a dictionary of :class:`~sapphire.api.Station` objects
                        for each station.
        :return: combined detector and station offsets for given station,
                 relative to the reference station.

        """
        offset_stations = station_numbers + [sn for sn in list(offsets.keys())
                                             if sn not in station_numbers]

        offset_matrix = zeros((len(offset_stations), len(offset_stations)))
        error_matrix = zeros((len(offset_stations), len(offset_stations)))

        for i, sn in enumerate(offset_stations):
            for j, ref_sn in enumerate(offset_stations):
                try:
                    o, e = offsets[sn].station_timing_offset(ref_sn,
                                                             midnight_ts)
                except Exception:
                    o, e = NO_STATION_OFFSET
                else:
                    if isnan(o) and isnan(e):
                        o, e = NO_STATION_OFFSET
                offset_matrix[i, j] = -o
                offset_matrix[j, i] = o
                error_matrix[i, j] = e ** 2
                error_matrix[j, i] = e ** 2

        ref_sn, predecessors = self.determine_best_reference(error_matrix,
                                                             station_numbers)

        best_offsets = {}
        for sn in station_numbers:
            best_offset = self._reconstruct_best_offset(
                predecessors, sn, ref_sn, station_numbers, offset_matrix)
            best_offsets[sn] = self._calculate_offsets(offsets[sn],
                                                       midnight_ts,
                                                       best_offset)
        return best_offsets

    def determine_best_reference(self, error_matrix, station_numbers):
        paths, predecessors = shortest_path(error_matrix, method='FW',
                                            directed=False,
                                            return_predecessors=True)
        n = len(station_numbers)
        # Only consider station in coincidence for reference
        total_errors = paths[:n, :n].sum(axis=1)
        best_reference = station_numbers[total_errors.argmin()]

        return best_reference, predecessors

    def _reconstruct_best_offset(self, predecessors, sn, ref_sn,
                                 station_numbers, offset_matrix):
        offset = 0.
        if sn != ref_sn:
            i = station_numbers.index(sn)
            j = station_numbers.index(ref_sn)
            while j != i:
                prev_j = j
                j = predecessors[i, j]
                offset += offset_matrix[prev_j, j]
        return offset

    def _calculate_offsets(self, station, ts0, offset):
        """Calculate combined station and detector offsets

        :param station: :class:`~sapphire.api.Station` object.
        :param ts0: gps timestamp for which the offsets are valid.
        :param offset: station offset to a reference station.
        :return: combined detector and station offsets for given station,
                 relative to the reference station.

        """
        detector_offsets = station.detector_timing_offset(ts0)
        return [offset + d_off for d_off in detector_offsets]

    def __repr__(self):
        return ("<%s, cluster: %r, direct: %r, fit: %r, curved: %r>" %
                (self.__class__.__name__, self.cluster, self.direct, self.fit,
                 self.curved))


class CoincidenceDirectionReconstructionDetectors(
        CoincidenceDirectionReconstruction):

    """Reconstruct direction for coincidences using each detector

    Instead of only the first arrival time per station this class
    uses the arrival time in each detector for the reconstruction.

    """

    def reconstruct_coincidence(self, coincidence_events, station_numbers=None,
                                offsets=None, initial=None):
        """Reconstruct a single coincidence

        :param coincidence_events: a coincidence list consisting of one
                                   or more (station_number, event) tuples.
        :param station_numbers: list of station numbers, to only use
                                events from those stations.
        :param offsets: dictionary with detector offsets for each station.
                        These detector offsets should be relative to one
                        detector from a specific station.
        :param initial: dictionary with already fitted shower parameters.
        :return: list of theta, phi, and station numbers.

        """
        if len(coincidence_events) < 1:
            return nan, nan, []
        if offsets is None:
            offsets = {}
        if initial is None:
            initial = {}

        # Subtract base timestamp to prevent loss of precision
        ts0 = int(coincidence_events[0][1]['timestamp'])
        ets0 = ts0 * int(1e9)
        self.cluster.set_timestamp(ts0)
        t, x, y, z, nums = ([], [], [], [], [])

        offsets = self.get_station_offsets(coincidence_events, station_numbers,
                                           offsets, ts0)

        for station_number, event in coincidence_events:
            if station_numbers is not None:
                if station_number not in station_numbers:
                    continue
            t_off = offsets.get(station_number, NO_OFFSET)
            station = self.cluster.get_station(station_number)
            t_detectors = relative_detector_arrival_times(event, ets0,
                                                          offsets=t_off,
                                                          station=station)
            for t_detector, detector in zip(t_detectors, station.detectors):
                if not isnan(t_detector):
                    dx, dy, dz = detector.get_coordinates()
                    t.append(t_detector)
                    x.append(dx)
                    y.append(dy)
                    z.append(dz)
            if not all(isnan(t_detectors)):
                nums.append(station_number)

        if len(t) >= 3 and 'core_x' in initial and 'core_y' in initial:
            theta, phi = self.curved.reconstruct_common(t, x, y, z, initial)
        elif len(t) == 3:
            theta, phi = self.direct.reconstruct_common(t, x, y, z, initial)
        elif len(t) > 3:
            theta, phi = self.fit.reconstruct_common(t, x, y, z, initial)
        else:
            theta, phi = (nan, nan)

        return theta, phi, nums


class BaseDirectionAlgorithm(object):

    """No actual direction reconstruction algorithm

    Simply returns (nan, nan) as direction.

    """

    @classmethod
    def reconstruct_common(cls, t, x, y, z=None, initial=None):
        """Reconstruct shower angles

        :param t: detector arrival time in ns.
        :param x,y: positions of detectors in m.
        :param z: height of detectors in m.
        :param initial: dictionary containing values from previous
                        reconstructions.
        :return: reconstructed theta and phi angles.

        """
        return cls.reconstruct()

    @staticmethod
    def reconstruct():
        """Reconstruct shower angles

        :return: reconstructed theta and phi angles.

        """
        return (nan, nan)


class DirectAlgorithm(BaseDirectionAlgorithm):

    """Reconstruct angles using direct analytical formula.

    This implements the equations derived in Fokkema2012 sec 4.2.
    (DOI: 10.3990/1.9789036534383)

    This algorithm assumes each detector is at the same altitude.

    Note! The detectors are 0-based.

    """

    @classmethod
    def reconstruct_common(cls, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t: arrival times in detector 0, 1 and 2 in ns.
        :param x,y: positions of detector 0, 1 and 2 in m.
        :param z: height of detectors 0, 1 and 2 is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions is ignored.

        """
        if len(t) > 3 or len(x) > 3 or len(y) > 3:
            warning_only_three()

        dt = make_relative(t)
        dx = make_relative(x)
        dy = make_relative(y)

        r1 = vector_length(dx[1], dy[1])
        r2 = vector_length(dx[2], dy[2])

        phi1 = arctan2(dy[1], dx[1])
        phi2 = arctan2(dy[2], dx[2])

        return cls.reconstruct(dt[1], dt[2], r1, r2, phi1, phi2)

    @classmethod
    def reconstruct(cls, dt1, dt2, r1, r2, phi1, phi2):
        """Reconstruct angles from 3 detections

        :param dt#: arrival times in detector 1 and 2 relative to
                    detector 0 in ns (!).
        :param r#,phi#: position of detector 1 and 2 relative to
                        detector 0 in m and radians.
        :return: theta as given by Fokkema2012 eq 4.14,
                 phi as given by Fokkema2012 eq 4.13.

        """
        if dt1 == 0 and dt2 == 0:
            # No time difference means shower came from zenith.
            return 0, 0

        phi = arctan2(-(r1 * dt2 * cos(phi1) - r2 * dt1 * cos(phi2)),
                      (r1 * dt2 * sin(phi1) - r2 * dt1 * sin(phi2)))

        # The directional vector c * dt should be negative,
        # not apparent in Fokkema2012 fig 4.4.
        theta = nan
        if r1 == 0 or r2 == 0:
            pass
        elif not dt1 == 0 and not phi - phi1 == pi / 2:
            sintheta = c * -dt1 / (r1 * cos(phi - phi1))
            if abs(sintheta) <= 1:
                theta = arcsin(sintheta)
        elif not dt2 == 0 and not phi - phi2 == pi / 2:
            sintheta = c * -dt2 / (r2 * cos(phi - phi2))
            if abs(sintheta) <= 1:
                theta = arcsin(sintheta)

        # We limit theta to positive values.  If theta is negative, we
        # make it positive, but need to rotate phi by 180 degrees.
        if isnan(theta):
            phi = nan
        elif theta < 0:
            theta *= -1
            phi += pi
            phi = norm_angle(phi)

        return theta, phi

    @classmethod
    def rel_theta1_errorsq(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.23"""

        sintheta = sin(theta)
        sinphiphi1 = sin(phi - phi1)

        den = r1 ** 2 * (1 - sintheta ** 2) * cos(phi - phi1) ** 2

        aa = (r1 ** 2 * sinphiphi1 ** 2 *
              cls.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        bb = -(2 * r1 * c * sinphiphi1 *
               (cls.dphi_dt0(theta, phi, phi1, phi2, r1, r2) -
                cls.dphi_dt1(theta, phi, phi1, phi2, r1, r2)))
        cc = 2 * c ** 2

        errsq = (aa * sintheta ** 2 + bb * sintheta + cc) / den

        return where(isnan(errsq), inf, errsq)

    @classmethod
    def rel_theta2_errorsq(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.23"""

        sintheta = sin(theta)
        sinphiphi2 = sin(phi - phi2)

        den = r2 ** 2 * (1 - sintheta ** 2) * cos(phi - phi2) ** 2

        aa = (r2 ** 2 * sinphiphi2 ** 2 *
              cls.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        bb = -(2 * r2 * c * sinphiphi2 *
               (cls.dphi_dt0(theta, phi, phi1, phi2, r1, r2) -
                cls.dphi_dt2(theta, phi, phi1, phi2, r1, r2)))
        cc = 2 * c ** 2

        errsq = (aa * sintheta ** 2 + bb * sintheta + cc) / den

        return where(isnan(errsq), inf, errsq)

    @staticmethod
    def rel_phi_errorsq(theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.22"""

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) ** 2 * r1 ** 2 * r2 ** 2 * sin(theta) ** 2 *
               (sinphi1 * cos(phi - phi2) - sinphi2 * cos(phi - phi1)) ** 2 /
               c ** 2)

        aa = (r1 ** 2 * sinphi1 ** 2 +
              r2 ** 2 * sinphi2 ** 2 -
              r1 * r2 * sinphi1 * sinphi2)
        bb = (2 * r1 ** 2 * sinphi1 * cosphi1 +
              2 * r2 ** 2 * sinphi2 * cosphi2 -
              r1 * r2 * (sinphi2 * cosphi1 + sinphi1 * cosphi2))
        cc = (r1 ** 2 * cosphi1 ** 2 +
              r2 ** 2 * cosphi2 ** 2 -
              r1 * r2 * cosphi1 * cosphi2)

        return 2 * (aa * tanphi ** 2 + bb * tanphi + cc) / den

    @classmethod
    def dphi_dt0(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.19"""

        return -(cls.dphi_dt1(theta, phi, phi1, phi2, r1, r2) +
                 cls.dphi_dt2(theta, phi, phi1, phi2, r1, r2))

    @staticmethod
    def dphi_dt1(theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.20"""

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta) *
               (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2)) /
               c)
        num = -r2 * (sinphi2 * tanphi + cosphi2)

        return num / den

    @staticmethod
    def dphi_dt2(theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.21"""

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta) *
               (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2)) /
               c)
        num = r1 * (sinphi1 * tanphi + cosphi1)

        return num / den


class DirectAlgorithmCartesian(BaseDirectionAlgorithm):

    """Reconstruct angles using direct analytical formula.

    This implements the equations derived in Montanus2014.
    "Direction reconstruction of cosmic air showers with
    detectorstations at different altitudes"

    Here the 2D version is used, assuming each detector is at the same
    altitude.

    """

    @classmethod
    def reconstruct_common(cls, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t: arrival times in detector 0, 1 and 2 in ns.
        :param x,y: positions of detector 0, 1 and 2 in m.
        :param z: height of detectors 0, 1 and 2 is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions is ignored.

        """
        if len(t) > 3 or len(x) > 3 or len(y) > 3:
            warning_only_three()

        dt = make_relative(t)
        dx = make_relative(x)
        dy = make_relative(y)

        return cls.reconstruct(dt[1], dt[2], dx[1], dx[2], dy[1], dy[2])

    @staticmethod
    def reconstruct(dt1, dt2, dx1, dx2, dy1, dy2):
        """Reconstruct angles from 3 detections

        :param dt#: arrival times in detector 1 and 2 relative to
                    detector 0 in ns.
        :param dx#,dy#: position of detector 1 and 2 relative to
                         detector 0 in m.
        :return: theta as given by Montanus2014 eq 27,
                 phi as given by Montanus2014 eq 26.

        """
        ux = c * (dt2 * dx1 - dt1 * dx2)
        uy = c * (dt2 * dy1 - dt1 * dy2)

        vz = dx1 * dy2 - dx2 * dy1

        theta = nan
        phi = nan

        if not vz == 0:
            usquared = ux * ux + uy * uy
            vzsquared = vz * vz
            uvzsqrt = sqrt(usquared / vzsquared)
            if uvzsqrt <= 1.0:
                theta = arcsin(uvzsqrt)
                phi = arctan2(-ux * vz, uy * vz)

        return theta, phi


class DirectAlgorithmCartesian3D(BaseDirectionAlgorithm):

    """Reconstruct angles using direct analytical formula.

    This implements the equations derived in Montanus2014.
    "Direction reconstruction of cosmic air showers with
    detectorstations at different altitudes"

    Here the 3D version is used, assuming each detector is at the same
    altitude.

    """

    @classmethod
    def reconstruct_common(cls, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t: arrival times in detector 0, 1 and 2 in ns.
        :param x,y,z: positions of detector 0, 1 and 2 in m.
        :param initial: dictionary containing values from previous
                        reconstructions is ignored.

        """
        if z is None:
            z = [0] * len(x)

        if len(t) > 3 or len(x) > 3 or len(y) > 3 or len(z) > 3:
            warning_only_three()

        dt = make_relative(t)
        dx = make_relative(x)
        dy = make_relative(y)
        dz = make_relative(z)

        return cls.reconstruct(dt[1], dt[2], dx[1], dx[2], dy[1], dy[2], dz[1],
                               dz[2])

    @staticmethod
    def reconstruct(dt1, dt2, dx1, dx2, dy1, dy2, dz1=0, dz2=0):
        """Reconstruct angles from 3 detections

        :param dt#: arrival times in detector 1 and 2 relative to
                    detector 0 in ns.
        :param dx#,dy#,dz#: position of detector 1 and 2 relative to
                            detector 0 in m.
        :return: theta as given by Montanus2014 eq 24,
                 phi as given by Montanus2014 eq 22.

        """
        d1 = array([dx1, dy1, dz1])
        d2 = array([dx2, dy2, dz2])
        u = c * (dt2 * d1 - dt1 * d2)
        v = cross(d1, d2)
        uxv = cross(u, v)

        usquared = dot(u, u)
        vsquared = dot(v, v)
        underroot = vsquared - usquared

        theta = nan
        phi = nan

        if underroot > 0 and not vsquared == 0:
            term = v * sqrt(underroot)
            nplus = (uxv + term) / vsquared
            nmin = (uxv - term) / vsquared

            phiplus = arctan2(nplus[1], nplus[0])
            thetaplus = arccos(nplus[2])

            phimin = arctan2(nmin[1], nmin[0])
            thetamin = arccos(nmin[2])

            if isnan(thetaplus):
                thetaplus = pi

            if isnan(thetamin):
                thetamin = pi

            # Allow solution only if it is the only one above horizon
            if thetaplus <= pi / 2. and thetamin > pi / 2.:
                theta = thetaplus
                phi = phiplus
            elif thetaplus > pi / 2. and thetamin <= pi / 2.:
                theta = thetamin
                phi = phimin

        return theta, phi


class SphereAlgorithm(object):

    """Reconstruct the direction in equatorial coordinates

    Note: currently incompatible with the other algorithms!

    This class uses a different coordinate systems than the other
    algorithms. The location input is in ECEF coordinates and a
    timestamp is required to connect the direction to the equatorial
    coordinates.

    """

    @classmethod
    def reconstruct_equatorial(cls, t, x, y, z, timestamp):
        """Reconstructs the source in the Equatorial Coordinate System.

        :param t: An array with three arrival times in ns.
        :param x,y,z: arrays with the ECEF locations of the
                      three detectors / stations in meters.
        :param timestamp: The UTC timestamp of the coincidence in s.
        :return: declination and right ascension of the source. The
                 apparent location of the cosmic ray source in the
                 Equatorial Coordinate System.

        """
        t_int = array([-1000, -10000]) + t[0]
        x_int, y_int, z_int = cls.interaction_curve(x, y, z, t, t_int)
        dec = arctan2(z_int[1] - z_int[0],
                      sqrt((x_int[1] - x_int[0]) ** 2. +
                           (y_int[1] - y_int[0]) ** 2.))
        ra = arctan2(x_int[1] - x_int[0], y_int[1] - y_int[0])
        return dec, ra

    @staticmethod
    def interaction_curve(x, y, z, t, t_int):
        """Calculates the curve of possible primary interactions

        This uses the arrival times in three detectors. The algorithm is
        based on location calculations used for LORAN, DECCA, RACAL, GPS
        as described by N.G. Schultheiss 2012

        :param x,y,z: arrays with the orthogonal coordinates of the three
                      detection points in m.
        :param t: arrival times of the detectors in ns.
        :param t_int: interaction time in ns.
        :return: parameters x_int, y_int, z_int

        """
        x01 = x[0] - x[1]
        x02 = x[0] - x[2]
        y01 = y[0] - y[1]
        y02 = y[0] - y[2]
        z01 = z[0] - z[1]
        z02 = z[0] - z[2]
        t01 = t[0] - t[1]
        t02 = t[0] - t[2]

        a = 2. * (x01 * y02 - x02 * y01)
        b = 2. * (x02 * z01 - x01 * z02)
        h = 2. * (x02 * t01 - x01 * t02) * c ** 2
        d = (x02 * (x01 ** 2 + y01 ** 2 + z01 ** 2 - (t01 * c) ** 2) -
             x01 * (x02 ** 2 + y02 ** 2 + z02 ** 2 - (t02 * c) ** 2))
        e = 2. * (y01 * z02 - y02 * z01)
        f = 2. * (y01 * t02 - y02 * t01) * c ** 2
        g = (y01 * (x02 ** 2 + y02 ** 2 + z02 ** 2 - (t02 * c) ** 2) -
             y02 * (x01 ** 2 + y01 ** 2 + z01 ** 2 - (t01 * c) ** 2))

        t = a ** 2 + b ** 2 + e ** 2
        v = (b * h + e * f) / t
        w = (b * d + e * g) / t
        p = (d ** 2 + g ** 2) / t
        q = 2 * (h * d + f * g) / t
        r = (h ** 2 + f ** 2 - (a * c) ** 2) / t

        t_int0 = t_int - t[0]

        sign = 1

        z = -v * t_int0 - w + sign * sqrt((v ** 2 - r) * t_int0 ** 2 +
                                          (2 * v * w - q) * t_int0 +
                                          w ** 2 - p)
        y = (b * z + h * t_int0 + d) / a
        x = (e * z + f * t_int0 + g) / a

        x_int = x[0] + x
        y_int = y[0] + y
        z_int = z[0] + z

        int_length = vector_length(x_int[0], y_int[0], z_int[0])
        det_length = vector_length(x[0], y[0], z[0])

        if det_length > int_length:
            # Select interaction above the earths surface.

            sign = -1
            z = -v * t_int0 - w + sign * sqrt((v ** 2 - r) * t_int0 ** 2 +
                                              (2 * v * w - q) * t_int0 +
                                              w ** 2 - p)
            y = (b * z + h * t_int0 + d) / a
            x = (e * z + f * t_int0 + g) / a

            x_int = x[0] + x
            y_int = y[0] + y
            z_int = z[0] + z

        return x_int, y_int, z_int, t_int


class FitAlgorithm3D(BaseDirectionAlgorithm):

    @classmethod
    def reconstruct_common(cls, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 or more detections

        This function converts the arguments to be suitable for the
        algorithm.

        :param t: arrival times of the detectors in ns.
        :param x,y,z: positions of the detectors in m. The height
                      for all detectors will be set to 0 if not given.
        :param initial: dictionary containing values from previous
                        reconstructions is ignored.

        """
        if z is None:
            z = [0] * len(x)

        return cls.reconstruct(t, x, y, z)

    @classmethod
    def reconstruct(cls, t, x, y, z):
        """Reconstruct angles for many detections

        :param t#: arrival times in the detectors in ns.
        :param x#,y#,z#: position of the detectors in m.
        :return: theta as given by Montanus2014 eq 21,
                 phi as given by Montanus2014 eq 22.

        """
        if not logic_checks(t, x, y, z):
            return nan, nan

        dt = make_relative(t)[1:]
        dx = make_relative(x)[1:]
        dy = make_relative(y)[1:]
        dz = make_relative(z)[1:]

        cons = {'type': 'eq', 'fun': cls.constraint_normal_vector}

        fit = minimize(cls.best_fit, x0=(0.1, 0.1, 0.989, 0.),
                       args=(dt, dx, dy, dz), method="SLSQP",
                       bounds=((-1, 1), (-1, 1), (-1, 1), (None, None)),
                       constraints=cons,
                       options={'ftol': 1e-9, 'eps': 1e-7, 'maxiter': 50})
        if fit.success:
            phi1 = arctan2(fit.x[1], fit.x[0])
            theta1 = arccos(fit.x[2])
        else:
            phi1 = nan
            theta1 = nan

        fit = minimize(cls.best_fit, x0=(-0.1, -0.1, -0.989, 0.),
                       args=(dt, dx, dy, dz), method="SLSQP",
                       bounds=((-1, 1), (-1, 1), (-1, 1), (None, None)),
                       constraints=cons,
                       options={'ftol': 1e-9, 'eps': 1e-7, 'maxiter': 50})
        if fit.success:
            phi2 = arctan2(fit.x[1], fit.x[0])
            theta2 = arccos(fit.x[2])
        else:
            phi2 = nan
            theta2 = nan

        # In case one of the theta's is smaller than pi/2 (shower from above)
        # and the other is either nan or larger than pi/2 (shower from below),
        # the first one is considered correct.
        # If both come from above (or from below), both theta's are rejected.
        if theta1 <= pi / 2. and (isnan(theta2) or theta2 > pi / 2.):
            theta = theta1
            phi = phi1
        elif (isnan(theta1) or theta1 > pi / 2.) and theta2 <= pi / 2.:
            theta = theta2
            phi = phi2
        else:
            theta = nan
            phi = nan

        return theta, phi

    @staticmethod
    def constraint_normal_vector(n):
        """This should be equal to zero"""

        return n[0] ** 2 + n[1] ** 2 + n[2] ** 2 - 1

    @staticmethod
    def best_fit(n_xyz, dt, dx, dy, dz):
        """The function to be minimized to find the direction

        :param n_xyz: list containing the unit vector.
        :param dt: list of relative arrival times in the detectors in ns.
        :param dx,dy,dz: list of relative detector positions in m.
        :return: least sum of squares as in Montanus2014, eq 36

        """
        nx, ny, nz, m = n_xyz

        slq = sum([(nx * xi + ny * yi + zi * nz + c * ti + m) ** 2
                   for ti, xi, yi, zi in zip(dt, dx, dy, dz)])
        return slq + m * m


class RegressionAlgorithm(BaseDirectionAlgorithm):

    """Reconstruct angles using an analytical regression formula.

    This implements the equations as for ISVHECRI (Montanus 2014).
    "Direction reconstruction of cosmic air showers with
    three or more detectorstations in a horizontal (for the
    moment) plane"

    """

    @classmethod
    def reconstruct_common(cls, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 or more detections

        This function converts the arguments to be suitable for the
        algorithm.

        :param t: arrival times of the detectors in ns.
        :param x,y,z: positions of the detectors in m. The height
                      is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions is ignored.

        """
        return cls.reconstruct(t, x, y)

    @classmethod
    def reconstruct(cls, t, x, y):
        """Reconstruct angles for many detections

        :param t: arrival times in the detectors in ns.
        :param x,y: positions of the detectors in m.
        :return: theta as derived by Montanus2014,
                 phi as derived by Montanus2014.

        """
        if not logic_checks(t, x, y, [0] * len(t)):
            return nan, nan

        k = len(t)
        xs = sum(x)
        ys = sum(y)
        ts = sum(t)

        xx = 0.
        yy = 0.
        tx = 0.
        ty = 0.
        xy = 0.

        for ti, xi, yi in zip(t, x, y):
            xx += xi ** 2
            yy += yi ** 2
            tx += ti * xi
            ty += ti * yi
            xy += xi * yi

        denom = (k * xy ** 2 + xs ** 2 * yy + ys ** 2 * xx - k * xx * yy -
                 2 * xs * ys * xy)
        if denom == 0:
            denom = nan

        numer = (tx * (k * yy - ys ** 2) + xy * (ts * ys - k * ty) +
                 xs * ys * ty - ts * xs * yy)
        nx = c * numer / denom

        numer = (ty * (k * xx - xs ** 2) + xy * (ts * xs - k * tx) +
                 xs * ys * tx - ts * ys * xx)
        ny = c * numer / denom

        horiz = nx ** 2 + ny ** 2
        if horiz > 1.:
            theta = nan
            phi = nan
        else:
            nz = sqrt(1 - nx ** 2 - ny ** 2)
            phi = arctan2(ny, nx)
            theta = arccos(nz)

        return theta, phi


class RegressionAlgorithm3D(BaseDirectionAlgorithm):

    """Reconstruct angles by iteratively applying a regression formula.

    This implements the equations as recently derived (Montanus 2014).
    "Direction reconstruction of cosmic air showers with
    three or more detectorstations at arbitrary altitudes"

    """

    MAX_ITERATIONS = 1000

    @classmethod
    def reconstruct_common(cls, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 or more detections

        This function converts the arguments to be suitable for the
        algorithm.

        :param t: arrival times of the detectors in ns.
        :param x,y,z: positions of the detectors in m. The height
                      for all detectors will be set to 0 if not given.
        :param initial: dictionary containing values from previous
                        reconstructions is ignored.

        """
        if z is None:
            z = [0] * len(x)

        return cls.reconstruct(t, x, y, z)

    @classmethod
    def reconstruct(cls, t, x, y, z):
        """Reconstruct angles for many detections

        :param t: arrival times in the detectors in ns.
        :param x,y,z: positions of the detectors in m.
        :return: theta as derived by Montanus2014,
                 phi as derived by Montanus2014.

        """
        if not logic_checks(t, x, y, z):
            return nan, nan

        regress2d = RegressionAlgorithm()
        theta, phi = regress2d.reconstruct_common(t, x, y)

        dtheta = 1.
        iteration = 0
        while dtheta > 0.001:
            iteration += 1
            if iteration > cls.MAX_ITERATIONS:
                return nan, nan
            nxnz = tan(theta) * cos(phi)
            nynz = tan(theta) * sin(phi)
            nz = cos(theta)
            x_proj = [xi - zi * nxnz for xi, zi in zip(x, z)]
            y_proj = [yi - zi * nynz for yi, zi in zip(y, z)]
            t_proj = [ti + zi / (c * nz) for ti, zi in zip(t, z)]
            theta_prev = theta
            theta, phi = regress2d.reconstruct_common(t_proj, x_proj, y_proj)
            dtheta = abs(theta - theta_prev)

        return theta, phi


class CurvedMixin(object):

    """Provide methods to estimate the time delay due to front curvature

    Given a core location, detector position, and shower angle the radial core
    distance can be determined, which can be used to determine the expected
    time delay.

    """

    def time_delay(self, x, y, core_x, core_y, theta, phi):
        r = self.radial_core_distance(x, y, core_x, core_y, theta, phi)
        return self.front.delay_at_r(r)

    @classmethod
    def radial_core_distance(cls, x, y, core_x, core_y, theta, phi):
        """Determine the radial core distance

        :param x,y,z: positions of the detectors in m.
        :param core_x,core_y: core position at z = 0 in m.
        :param theta,phi: reconstructed shower direction.
        :return: radial core distance in m.

        """
        dx = core_x - x
        dy = core_y - y
        nx = sin(theta) * cos(phi)
        ny = sin(theta) * sin(phi)
        return sqrt(dx ** 2 * (1 - nx ** 2) + dy ** 2 * (1 - ny ** 2) -
                    2 * dx * dy * nx * ny)


class CurvedRegressionAlgorithm(CurvedMixin, BaseDirectionAlgorithm):

    """Reconstruct angles taking the shower front curvature into account.

    Take the shower front curvature into account. Assumes knowledge about the
    shower core position.

    """

    MAX_ITERATIONS = 1000

    def __init__(self):
        self.front = CorsikaStationFront()

    def reconstruct_common(self, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 or more detections

        This function converts the arguments to be suitable for the
        algorithm.

        :param t: arrival times of the detectors in ns.
        :param x,y,z: positions of the detectors in m.  The height
                      is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions, including core position.

        """
        if initial is None:
            initial = {}
        core_x = initial.get('core_x', nan)
        core_y = initial.get('core_y', nan)
        if isnan(core_y) or isnan(core_y):
            return nan, nan

        return self.reconstruct(t, x, y, core_x, core_y)

    def reconstruct(self, t, x, y, core_x, core_y):
        """Reconstruct angles for many detections

        :param t: arrival times in the detectors in ns.
        :param x,y: positions of the detectors in m.
        :param core_x,core_y: core position at z = 0 in m.
        :return: theta as derived by Montanus2014,
                 phi as derived by Montanus2014.

        """
        if not logic_checks(t, x, y, [0] * len(t)):
            return nan, nan

        regress2d = RegressionAlgorithm()
        theta, phi = regress2d.reconstruct_common(t, x, y)

        dtheta = 1.
        iteration = 0
        while dtheta > 0.001:
            iteration += 1
            if iteration > self.MAX_ITERATIONS:
                return nan, nan
            t_proj = [ti - self.time_delay(xi, yi, core_x, core_y, theta, phi)
                      for ti, xi, yi in zip(t, x, y)]
            theta_prev = theta
            theta, phi = regress2d.reconstruct_common(t_proj, x, y)
            dtheta = abs(theta - theta_prev)

        return theta, phi


class CurvedRegressionAlgorithm3D(CurvedMixin, BaseDirectionAlgorithm):

    """Reconstruct angles accounting for front curvature and detector altitudes

    Take the shower front curvature and different detector heights into
    account. Assumes knowledge about the shower core position.

    """

    MAX_ITERATIONS = 1000

    def __init__(self):
        self.front = CorsikaStationFront()

    def reconstruct_common(self, t, x, y, z=None, initial=None):
        """Reconstruct angles from 3 or more detections

        This function converts the arguments to be suitable for the
        algorithm.

        :param t: arrival times of the detectors in ns.
        :param x,y,z: positions of the detectors in m. The height
                      for all detectors will be set to 0 if not given.
        :param initial: dictionary containing values from previous
                        reconstructions, including core position.

        """
        if initial is None:
            initial = {}
        core_x = initial.get('core_x', nan)
        core_y = initial.get('core_y', nan)
        if isnan(core_y) or isnan(core_y):
            return nan, nan

        if z is None:
            z = [0] * len(x)

        return self.reconstruct(t, x, y, z, core_x, core_y)

    def reconstruct(self, t, x, y, z, core_x, core_y):
        """Reconstruct angles for many detections

        :param t: arrival times in the detectors in ns.
        :param x,y,z: positions of the detectors in m.
        :param core_x,core_y: core position at z = 0 in m.
        :return: theta as derived by Montanus2014,
                 phi as derived by Montanus2014.

        """
        if not logic_checks(t, x, y, z):
            return nan, nan

        regress2d = RegressionAlgorithm()
        theta, phi = regress2d.reconstruct_common(t, x, y)

        dtheta = 1.
        iteration = 0
        while dtheta > 0.001:
            iteration += 1
            if iteration > self.MAX_ITERATIONS:
                return nan, nan
            nxnz = tan(theta) * cos(phi)
            nynz = tan(theta) * sin(phi)
            nz = cos(theta)
            x_proj = [xi - zi * nxnz for xi, zi in zip(x, z)]
            y_proj = [yi - zi * nynz for yi, zi in zip(y, z)]
            t_proj = [ti + zi / (c * nz) -
                      self.time_delay(xpi, ypi, core_x, core_y, theta, phi)
                      for ti, xpi, ypi, zi in zip(t, x_proj, y_proj, z)]
            theta_prev = theta
            theta, phi = regress2d.reconstruct_common(t_proj, x_proj, y_proj)
            dtheta = abs(theta - theta_prev)

        return theta, phi


def logic_checks(t, x, y, z):
    """Check for impossible reconstructions

    Criteria:

    - No two detectors are at the same position.
    - Time difference between two detections should be less than distance / c.
    - All detectors on a line is bad.

    To fix:

    - Time difference can still be to large in cases where a different
      distance becomes relevant.

    :param t: arrival times in the detectors in ns.
    :param x,y,z: positions of the detectors in m.
    :return: True if the checks pass, False otherwise.

    """
    # Check for identical positions
    if len(t) == 3:
        xyz = list(zip(x, y, z))
        if not len(xyz) == len(set(xyz)):
            return False

    txyz = list(zip(t, x, y, z))

    # Check if the time difference it larger than expected by c
    if len(t) == 3:
        for txyz0, txyz1 in combinations(txyz, 2):
            dt = abs(txyz0[0] - txyz1[0])
            dx = txyz0[1] - txyz1[1]
            dy = txyz0[2] - txyz1[2]
            dz = txyz0[3] - txyz1[3]
            dt_max = vector_length(dx, dy, dz) / c
            if dt_max < dt:
                return False

    # Check if all the positions are (almost) on a single line
    largest_of_smallest_angles = 0
    for txyz0, txyz1, txyz2 in combinations(txyz, 3):
        dx1 = txyz0[1] - txyz1[1]
        dy1 = txyz0[2] - txyz1[2]
        dz1 = txyz0[3] - txyz1[3]
        dx2 = txyz0[1] - txyz2[1]
        dy2 = txyz0[2] - txyz2[2]
        dz2 = txyz0[3] - txyz2[3]
        dx3 = dx2 - dx1
        dy3 = dy2 - dy1
        dz3 = dz2 - dz1
        lenvec01 = vector_length(dx1, dy1, dz1)
        lenvec02 = vector_length(dx2, dy2, dz2)
        lenvec12 = vector_length(dx3, dy3, dz3)

        # area triangle is |cross product|
        area = abs(dx1 * dy2 - dx2 * dy1 + dy1 * dz2 - dy2 * dz1 +
                   dz1 * dx2 - dz2 * dx1)

        # prevent floating point errors
        if area < 1e-7:
            return False

        # sine of angle is area divided by two sides
        sin1 = area / lenvec01 / lenvec02
        sin2 = area / lenvec01 / lenvec12
        sin3 = area / lenvec02 / lenvec12

        # smallest sine
        smallest_angle = min(sin1, sin2, sin3)

        # remember largest of smallest sines
        largest_of_smallest_angles = max(largest_of_smallest_angles,
                                         smallest_angle)

    # discard reconstruction if the largest of the smallest angles of each
    # triangle is smaller than 0.1 rad (5.73 degrees)
    if largest_of_smallest_angles < 0.1:
        return False

    return True


def warning_only_three():
    warnings.warn('Only the first three detections will be used')
