"""Perform simple simulations for timing

Throw a flat shower front on a cluster from various angles.
Simulate just the arrival times.

Example usage::

    import tables

    from sapphire.simulations.showerfront import FlatFrontSimulation
    from sapphire.clusters import ScienceParkCluster

    data = tables.open_file('/tmp/test_showerfront_simulation.h5', 'w')
    cluster = ScienceParkCluster()

    sim = FlatFrontSimulation(cluster, data, '/', 200)
    sim.run()

"""
from math import pi, sin, cos, sqrt, acos, floor

import numpy as np

from .base import BaseSimulation
from ..utils import pbar


class FlatFrontSimulation(BaseSimulation):

    def __init__(self, *args, **kwargs):
        super(FlatFrontSimulation, self).__init__(*args, **kwargs)

        for station in self.cluster.stations:
            station.gps_offset = np.random.normal(0, 16)
            # The actual distribution is not yet very clear.
            # We assume it is gaussian for convenience.
            # Then the stddev is about 16 ns.

        # Store updated version of the cluster
        self.coincidence_group._v_attrs.cluster = self.cluster

    def __del__(self):
        self.finish()

    def finish(self):
        """Clean-up after simulation"""

        self.corsikafile.close()

    def generate_shower_parameters(self):
        """Generate shower parameters, i.e. azimuth and zenith angles.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :returns: dictionary with shower parameters: core_pos
                  (x, y-tuple) and azimuth.

        """
        shower_parameters = {'core_pos': (None, None),
                             'zenith': None,
                             'azimuth': None,
                             'size': None,
                             'energy': None,
                             'ext_timestamp': None}

        # DF: This is the fastest implementation I could come up with.  I
        # timed several permutations of numpy / math, and tried a Monte
        # Carlo method in which I pick x and y in a square (fast) and then
        # determine if they fall inside a circle or not (surprisingly
        # slow, because of an if-statement, and despite some optimizations
        # suggested by HM).

        for i in pbar(range(self.N)):
            azimuth = np.random.uniform(-pi, pi)
            zenith = self.inverse_zenith_probability(np.random.random())

            shower_parameters = {'ext_timestamp': (long(1e9) + i) * long(1e9),
                                 'azimuth': azimuth,
                                 'zenith': zenith,
                                 'core_pos': (None, None),
                                 'size': None,
                                 'energy': None}
            self._prepare_cluster_for_shower(azimuth)

            yield shower_parameters

    def inverse_zenith_probability(self, p):
        """Inverse cumulative probability distribution for zenith

        Derrived from Schultheiss "The acceptancy of the HiSPARC Network",
        (internal note), eq 2.4 from Rossi.

        :param p: probability value between 0 and 1.
        :returns: zenith with corresponding cumulative probability.

        """
        return acos((1 - p) ** (1 / 8.))

    def _prepare_cluster_for_shower(self, alpha):
        """Prepare the cluster object for the simulation of a shower.

        Rotate the cluster so that 'East' coincides with the shower
        azimuth direction.

        """
        # rotate the around the cluster center
        self.cluster.set_coordinates(0, 0, 0, -alpha)

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Return the arrival time of shower front passing the center of
        the detector.

        """
        arrival_time = (self.get_arrival_time(detector, shower_parameters) +
                        self.simulate_signal_transport_time(1)[0])
        observables = {'t': arrival_time}

        return observables

    def get_arrival_time(self, detector, shower_parameters):
        """Calculate arrival time

        Assumes a flat shower front and core position to be
        the center of the cluster. Does not take detector
        altitude into account.

        Equation from Fokkema2012 sec 4.2, eq 4.9.
        (DOI: 10.3990/1.9789036534383)

        :returns: Shower front arrival time in ns.

        """
        c = 0.3
        r1, phi1 = detector.get_polar_coordinates()
        phi = shower_parameters['azimuth']
        theta = shower_parameters['zenith']
        cdt = r1 * cos(phi - phi1) * sin(theta)
        return cdt / c

    def simulate_signal_transport_time(self, size):
        """ Simulate transport times of scintillation light to the PMT

        Generates random transit times within a given distribution and adds it
        to the times the particles passed the detector.

        """
        numbers = np.random.random(size)
        dt = []

        for x in numbers:
            if  x < 0.39377:
                dt.append(2.5507 + 2.39885 * x)
            else:
                dt.append(1.56764 + 4.89536 * x)

        return dt

    def simulate_gps(self, station_observables, shower_parameters, station):
        """Simulate gps timestamp

        Ensure that all detector arrival times are positive.

        """
        n_detectors = len(station.detectors)
        ids = range(1, n_detectors + 1)
        arrival_times = [station_observables['t%d' % id] for id in ids]
        ext_timestamp = shower_parameters['ext_timestamp']

        first_time = floor(sorted(arrival_times)[0])
        for id in ids:
            station_observables['t%d' % id] -= first_time

        arrival_times = [station_observables['t%d' % id] for id in ids]
        trigger_time = sorted(arrival_times)[1]
        station_observables['t_trigger'] = trigger_time

        ext_timestamp += int(first_time + trigger_time + station.gps_offset +
                             self.simulate_gps_uncertainty())
        timestamp = int(ext_timestamp / long(1e9))
        nanoseconds = int(ext_timestamp % long(1e9))

        gps_timestamp = {'ext_timestamp': ext_timestamp,
                         'timestamp': timestamp, 'nanoseconds': nanoseconds}
        station_observables.update(gps_timestamp)

        return station_observables

    def simulate_gps_uncertainty(self):
        """Simulate uncertainty from GPS receiver"""

        return np.random.normal(0, 4.5)
