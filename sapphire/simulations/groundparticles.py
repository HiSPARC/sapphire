from math import pi, sin, cos, sqrt
import warnings

import numpy as np
import tables
import progressbar

from .base import BaseSimulation


class GroundParticlesSimulation(BaseSimulation):

    def __init__(self, corsikafile_path, max_core_distance, *args, **kwargs):
        super(GroundParticlesSimulation, self).__init__(*args, **kwargs)

        self.corsikafile = tables.openFile(corsikafile_path, 'r')
        self.groundparticles = self.corsikafile.getNode('/groundparticles')
        self.max_core_distance = max_core_distance

        for station in self.cluster.stations:
            station.gps_offset = np.random.normal(0, 18)
            # The distribution and spread has to be determined more accurate.
            # At a later stage this has to be adjusted.

        self.coincidence_group._v_attrs.cluster = self.cluster

    def finish(self):
        """Clean-up after simulation"""

        self.corsikafile.close()

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :returns: dictionary with shower parameters: core_pos
                  (x, y-tuple) and azimuth.

        """
        R = self.max_core_distance

        event_header = self.groundparticles._v_attrs['event_header']
        corsika_parameters = {'zenith': event_header.zenith,
                              'energy': event_header.energy,
                              'particle': event_header.particle}

        pbar = progressbar.ProgressBar(widgets=[progressbar.Percentage(),
                                                progressbar.Bar(),
                                                progressbar.ETA()])

        # DF: This is the fastest implementation I could come up with.  I
        # timed several permutations of numpy / math, and tried a Monte
        # Carlo method in which I pick x and y in a square (fast) and then
        # determine if they fall inside a circle or not (surprisingly
        # slow, because of an if-statement, and despite some optimizations
        # suggested by HM).

        for i in pbar(range(self.N)):
            r = sqrt(np.random.uniform(0, R ** 2))
            phi = np.random.uniform(-pi, pi)
            alpha = np.random.uniform(-pi, pi)

            x = r * cos(phi)
            y = r * sin(phi)

            # Add Corsika shower azimuth to generated alpha
            # and make them fit in (-pi, pi].
            shower_azimuth = alpha + event_header.azimuth
            if shower_azimuth > pi:
                shower_azimuth -= 2 * pi
            elif shower_azimuth <= -pi:
                shower_azimuth += 2 * pi

            shower_parameters = {'ext_timestamp': (long(1e9) + i) * long(1e9),
                                 'core_pos': (x, y),
                                 'azimuth': shower_azimuth}
            shower_parameters.update(corsika_parameters)
            self._prepare_cluster_for_shower(x, y, alpha)

            yield shower_parameters

    def _prepare_cluster_for_shower(self, x, y, alpha):
        """Prepare the cluster object for the simulation of a shower.

        Rotate and translate the cluster so that (0, 0) coincides with the
        shower core position and 'East' coincides with the shower azimuth
        direction.

        """
        # rotate the core position around the original cluster center
        xp = x * cos(-alpha) - y * sin(-alpha)
        yp = x * sin(-alpha) + y * cos(-alpha)

        self.cluster.set_xyalpha_coordinates(-xp, -yp, -alpha)

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Checks if leptons have passed a detector. If so, it returns the number
        of leptons in the detector and the arrival time of the first lepton
        passing the detector.

        The detector is approximated by a square with a surface of 0.5
        square meter which is *not* correctly rotated.  In fact, during
        the simulation, the rotation of the detector is undefined.  This
        might be faster than a more thorough implementation.

        """
        detector_boundary = 0.3535534
        x, y = detector.get_xy_coordinates()

        # particle ids 2, 3, 5, 6 are electrons and muons, and id 4 is no
        # longer used (were neutrino's).
        query = ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
                 ' & (particle_id >= 2) & (particle_id <= 6)' %
                 (x - detector_boundary, x + detector_boundary,
                  y - detector_boundary, x + detector_boundary))
        detected = [row['t'] for row in self.groundparticles.where(query)]

        if detected:
            n_detected = len(detected)
            transporttimes = self.simulate_signal_transport_time(n_detected)
            for i in range(n_detected):
                detected[i] += transporttimes[i]
            observables = {'n': n_detected, 't': min(detected)}
        else:
            observables = {'n': 0., 't': -999}

        return observables

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

    def simulate_trigger(self, detector_observables):
        """Simulate a trigger response.

        :returns: True if at least 2 detectors detect at least one particle,
                  False otherwise.

        """
        detectors_hit = [True for observables in detector_observables
                         if observables['n'] > 0]

        if sum(detectors_hit) >= 2:
            return True
        else:
            return False

    def simulate_gps(self, station_observables, shower_parameters, station):
        """Simulate gps timestamp."""

        arrival_times = [station_observables['t%d' % id]
                         for id in range(1, 5)
                         if station_observables.get('n%d' % id, -1) > 0]

        if len(arrival_times) > 1:
            trigger_time = sorted(arrival_times)[1]

            timestamp = int(shower_parameters['ext_timestamp'] / 1000000000)
            nanoseconds = int(trigger_time + self.simulate_gps_uncertainty() +
                              station.gps_offset)
            ext_timestamp = shower_parameters['ext_timestamp'] + nanoseconds

            gps_timestamp = {'ext_timestamp': ext_timestamp,
                             'timestamp': timestamp, 'nanoseconds': nanoseconds}
            station_observables.update(gps_timestamp)

        return station_observables

    def simulate_gps_uncertainty(self):
        """Simulate uncertainty from GPS receiver"""

        return np.random.normal(0, 4.5)
