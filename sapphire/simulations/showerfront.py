"""Perform simple simulations for timing

Throw a simulated shower front on a cluster from various angles.
Simulate just the arrival times.

Example usage::

    >>> import tables

    >>> from sapphire import FlatFrontSimulation, ScienceParkCluster

    >>> data = tables.open_file('/tmp/test_showerfront_simulation.h5', 'w')
    >>> cluster = ScienceParkCluster()

    >>> sim = FlatFrontSimulation(cluster, data, '/', 200)
    >>> sim.run()

"""
from __future__ import division

from math import sin, cos, tan, atan2, sqrt

import numpy as np

from .detector import HiSPARCSimulation, ErrorlessSimulation
from ..utils import pbar, c, vector_length


class FlatFrontSimulation(HiSPARCSimulation):

    def __init__(self, *args, **kwargs):
        super(FlatFrontSimulation, self).__init__(*args, **kwargs)

        # Since the cluster is not rotated detector positions can be cached.
        for station in self.cluster.stations:
            for detector in station.detectors:
                detector.cylindrical_coordinates = \
                    detector.get_cylindrical_coordinates()

    def generate_shower_parameters(self):
        """Generate shower parameters, i.e. azimuth and zenith angles.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple) and azimuth.

        """
        for i in pbar(range(self.n), show=self.progress):
            shower_parameters = {'ext_timestamp': (int(1e9) + i) * int(1e9),
                                 'azimuth': self.generate_azimuth(),
                                 'zenith': self.generate_attenuated_zenith(),
                                 'core_pos': (None, None),
                                 'size': None,
                                 'energy': None}

            yield shower_parameters

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Return the arrival time of shower front passing the center of
        the detector.

        """
        arrival_time = self.simulate_adc_sampling(
            self.get_arrival_time(detector, shower_parameters) +
            self.simulate_signal_transport_time(1)[0] +
            detector.offset)
        observables = {'t': arrival_time}

        return observables

    def get_arrival_time(self, detector, shower_parameters):
        """Calculate arrival time

        Assumes a flat shower front and core position to be
        the center of the cluster.

        Equation based on Fokkema2012 sec 4.2, eq 4.9.
        With additions to account for altitude.
        (DOI: 10.3990/1.9789036534383)
        The directional vector c * dt should be negative,
        not apparent in Fokkema2012 fig 4.4.

        :return: Shower front arrival time in ns.

        """
        r1, phi1, z1 = detector.cylindrical_coordinates
        phi = shower_parameters['azimuth']
        theta = shower_parameters['zenith']
        r = r1 * cos(phi - phi1) - z1 * tan(theta)
        cdt = - (r * sin(theta) + z1 / cos(theta))
        dt = cdt / c
        return dt

    def simulate_gps(self, station_observables, shower_parameters, station):
        """Simulate gps timestamp

        Ensure that all detector arrival times are positive.

        """
        n_detectors = len(station.detectors)
        ids = list(range(1, n_detectors + 1))
        arrival_times = [station_observables['t%d' % id] for id in ids]
        ext_timestamp = shower_parameters['ext_timestamp']

        first_time = sorted(arrival_times)[0]
        for id in ids:
            station_observables['t%d' % id] -= first_time

        arrival_times = [station_observables['t%d' % id] for id in ids]
        trigger_time = sorted(arrival_times)[1]

        ext_timestamp += int(first_time + trigger_time + station.gps_offset +
                             self.simulate_gps_uncertainty())
        timestamp = int(ext_timestamp / int(1e9))
        nanoseconds = int(ext_timestamp % int(1e9))

        gps_timestamp = {'ext_timestamp': ext_timestamp,
                         'timestamp': timestamp,
                         'nanoseconds': nanoseconds,
                         't_trigger': trigger_time}
        station_observables.update(gps_timestamp)

        return station_observables


class FlatFrontSimulationWithoutErrors(ErrorlessSimulation,
                                       FlatFrontSimulation):

    """This simulation does not simulate errors/uncertainties

    This should result in perfect timing for the detectors.

    """

    pass


class FlatFrontSimulation2D(FlatFrontSimulation):

    """This simulation ignores detector altitudes."""

    def get_arrival_time(self, detector, shower_parameters):
        """Calculate arrival time

        Ignore detector altitudes

        Equation based on Fokkema2012 sec 4.2, eq 4.9.
        (DOI: 10.3990/1.9789036534383)

        """
        r1, phi1, _ = detector.cylindrical_coordinates
        phi = shower_parameters['azimuth']
        theta = shower_parameters['zenith']
        r = r1 * cos(phi - phi1)
        cdt = -r * sin(theta)
        dt = cdt / c
        return dt


class FlatFrontSimulation2DWithoutErrors(FlatFrontSimulation2D,
                                         FlatFrontSimulationWithoutErrors):

    """Ignore altitude of detectors and do not simulate errors."""

    pass


class ConeFrontSimulation(FlatFrontSimulation):

    """This simulation uses a cone shaped shower front.

    The opening angle of the cone is given in the init

    Example usage::

        >>> import tables

        >>> from sapphire import ConeFrontSimulation, ScienceParkCluster

        >>> data = tables.open_file('/tmp/test_showerfront_simulation.h5', 'w')
        >>> cluster = ScienceParkCluster()

        >>> sim = ConeFrontSimulation(100, cluster, data, '/', 200)
        >>> sim.run()

    """

    def __init__(self, max_core_distance, *args, **kwargs):
        """Calculate arrival time

        :param cone_angle: half of the opening angle of the cone.

        """
        super(ConeFrontSimulation, self).__init__(*args, **kwargs)
        self.max_core_distance = max_core_distance
        self.front = ConeFront()

    def generate_shower_parameters(self):
        """Generate shower parameters

        For this cone-shaped showerfront, the core position, the azimuth
        and zenith angle of the shower are generated.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple), azimuth and zenith.

        """
        r_max = self.max_core_distance

        for i in pbar(range(self.n), show=self.progress):
            x, y = self.generate_core_position(r_max)
            azimuth = self.generate_azimuth()

            shower_parameters = {'ext_timestamp': (int(1e9) + i) * int(1e9),
                                 'azimuth': azimuth,
                                 'zenith': self.generate_attenuated_zenith(),
                                 'core_pos': (x, y),
                                 'size': None,
                                 'energy': self.generate_energy(1e15, 1e17)}

            self._prepare_cluster_for_shower(x, y, azimuth)

            yield shower_parameters

    def _prepare_cluster_for_shower(self, x, y, alpha):
        """Prepare the cluster object for the simulation of a shower.

        Rotate and translate the cluster so that (0, 0, 0) coincides with the
        shower core position and 'East' coincides with the shower azimuth
        direction.

        """
        # rotate the core position around the original cluster center
        xp = x * cos(-alpha) - y * sin(-alpha)
        yp = x * sin(-alpha) + y * cos(-alpha)

        self.cluster.set_coordinates(-xp, -yp, 0, -alpha)

    def get_arrival_time(self, detector, shower_parameters):
        """Calculate arrival time"""

        x, y, z = detector.get_coordinates()
        r1 = vector_length(x, y)
        phi1 = atan2(y, x)

        phi = shower_parameters['azimuth']
        theta = shower_parameters['zenith']
        r = r1 * cos(phi - phi1) - z * tan(theta)
        cdt = - (r * sin(theta) + z / cos(theta))

        nx = sin(theta) * cos(phi)
        ny = sin(theta) * sin(phi)
        nz = cos(theta)

        r_core = sqrt(x ** 2 + y ** 2 + z ** 2 -
                      (x * nx + y * ny + z * nz) ** 2)
        t_shape = self.delay_at_r(r_core)
        dt = t_shape + (cdt / c)

        return dt


class FlatFront(object):

    """Simple flat shower front"""

    def delay_at_r(self, r):
        return 0.

    def front_shape(self, r):
        return 0.


class ConeFront(object):

    """Simple cone shaped shower front"""

    def delay_at_r(self, r):
        return self.front_shape(r)

    def front_shape(self, r):
        """Delay of the showerfront relative to flat as function of distance

        :param r: distance to the shower core in shower frame.
        :return: delay time of shower front.

        """
        return r * 0.2


class CorsikaStationFront(object):

    """Shower front shape derrived from CORSIKA simulations on a station.

    A set of CORSIKA generated showers were used to determine the median
    detected arrival time in a 4-detector station as a function of core
    distance.

    At large core distances the detection probability decreases and the
    arrival time becomes less accurate.

    Currently only support for energies between 1e15 and 1e17 eV.

    """

    def delay_at_r(self, r, energy=1e16, particle='proton'):
        return self.front_shape(r, energy, particle)

    def front_shape(self, r, energy, particle='proton'):
        if particle == 'proton':
            energies = [15, 15.5, 16, 16.5, 17]
            param_a = [0.0995, 0.0190, 0.0065, 0.0032, 0.00044]
            param_b = [1.05, 1.38, 1.58, 1.69, 2.01]
        elif particle == 'gamma':
            energies = [15, 15.5, 16, 17]
            param_a = [0.08466, 0.01277, 0.00433, 0.00067]
            param_b = [1.13, 1.51, 1.70, 2.00]
        elif particle == 'iron':
            energies = [15, 16, 17]
            param_a = [0.44666, 0.01301, 0.00039]
            param_b = [0.69, 1.39, 1.98]

        a = np.interp(np.log10(energy), energies, param_a)
        b = np.interp(np.log10(energy), energies, param_b)

        return self._front_shape(r, a, b)

    def _front_shape(self, r, a, b):
        return a * r ** b
