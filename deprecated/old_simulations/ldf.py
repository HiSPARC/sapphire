from scipy.special import gamma
import progressbar as pb
import sys

import numpy as np
from numpy import pi, arccos, sin, cos, sqrt

from base import BaseSimulation


class BaseLdfSimulation(BaseSimulation):
    """ Simulation using a lateral distribution function as a model

        Not using an EAS simulation but rather a lateral distribution function
        as the model for the detector simulation considerably speeds up the
        complete process.  However, we will have to validate the use of an LDF.

        This is the base class, so there is no LDF defined.  Override
        calculate_ldf_value to complete an implementation.

    """
    # Alternative for observables.nrows, which is only updated with a flush()
    # We don't want to flush often, so we use (and update!) this instead
    _observables_nrows = 0


    def __init__(self, cluster, data, output, R, N, shower_size=10 ** 4.8, **kwargs):
        self.shower_size = shower_size
        self.ldf = BaseLdf()

        super(BaseLdfSimulation, self).__init__(cluster, data, output, R, N, **kwargs)

    def run(self, positions=None, max_theta=None):
        """Run a simulation

        This is the code which performs the simulation.  It creates a list
        of positions, creates all necessary tables and performs the
        simulation.

        :param positions: if given, use these coordinates instead of
            generating new ones

        """
        self._run_welcome_msg()

        if not positions:
            positions = self.generate_positions()

        progress = pb.ProgressBar(maxval=self.N, widgets=[pb.Percentage(),
                                                          pb.Bar(),
                                                          pb.ETA()],
                                  fd=sys.stderr)

        for event_id, (r, phi) in progress(enumerate(positions)):
            shower_theta, shower_phi = self.generate_shower_theta_and_phi(max_theta)

            event = {'id': event_id, 'r': r, 'phi': phi, 'alpha': 0.,
                     'shower_theta': shower_theta, 'shower_phi': shower_phi,
                     'shower_size': self.shower_size}
            self.simulate_event(event)

        self.coincidences.flush()
        self.observables.flush()
        # integrity check
        assert(self._observables_nrows == self.observables.nrows)
        self.c_index.flush()

        self._run_exit_msg()

    def generate_shower_theta_and_phi(self, max_theta=None):
        if max_theta is None:
            shower_theta, shower_phi = 0., 0.
        else:
            min_x = cos(max_theta)
            x = np.random.uniform(min_x, 1.)
            shower_theta = arccos(x)
            shower_phi = np.random.uniform(-pi, pi)

        return shower_theta, shower_phi

    def simulate_event(self, event):
        multiplicity = 0
        station_event_ids = []

        for station in self.cluster.stations:
            has_trigger, id = self.simulate_station_observables_and_return_has_triggered_and_eventid(station, event)
            multiplicity += has_trigger
            station_event_ids.append(id)
        self.write_coincidence(event, multiplicity)
        self.c_index.append(station_event_ids)

    def simulate_station_observables_and_return_has_triggered_and_eventid(self, station, event):
        num_particles = []
        for detector in station.detectors:
            num_particles.append(self.simulate_detector_observables(detector, event))
        signals = self.simulate_detector_signals(num_particles)
        id = self.write_observables_and_return_id(station, event, signals)

        num_detectors_over_threshold = sum([True if u >= self.trig_threshold
                                            else False for u in signals])
        if num_detectors_over_threshold >= 2:
            has_triggered = True
        else:
            has_triggered = False

        return has_triggered, id

    def simulate_detector_observables(self, detector, event):
        R = self.calculate_core_distance(detector, event)
        density_on_front = self.calculate_ldf_value(R, event['shower_size'])
        density_on_ground = cos(event['shower_theta']) * density_on_front
        num_particles = density_on_ground * detector.get_area()

        return num_particles

    def calculate_core_distance(self, detector, event):
        r, phi = event['r'], event['phi']
        x0 = r * cos(phi)
        y0 = r * sin(phi)

        x1, y1 = detector.get_xy_coordinates()

        R = self.ldf.calculate_core_distance_from_coordinates_and_direction(
                x0, y0, x1, y1, event['shower_theta'], event['shower_phi'])

        return R

    def calculate_ldf_value(self, R, shower_size):
        return 0.

    def write_observables_and_return_id(self, station, event, num_particles):
        """Write observables from a single event

        :param station: Station instance
        :param n1, n2, n3, n4: number of particles detected for each detector
            1, 2, 3 and 4

        """
        row = self.observables.row

        x, y, alpha = station.get_xyalpha_coordinates()
        r, phi, alpha = station.get_rphialpha_coordinates()

        station = {'id': event['id'], 'station_id': station.station_id,
                   'r': r, 'phi': phi, 'alpha': alpha}

        timings = 4 * [0.]

        super(BaseLdfSimulation, self).write_observables(station, num_particles,
                                                         timings)

        self._observables_nrows += 1
        return self._observables_nrows - 1

    def write_coincidence(self, event, N):
        """Write coincidence information

        :param event: simulated shower event information
        :param N: number of stations which triggered

        """
        row = self.coincidences.row

        r = event['r']
        phi = event['phi']
        row['id'] = event['id']
        row['N'] = N
        row['r'] = r
        row['phi'] = phi
        row['x'] = r * cos(phi)
        row['y'] = r * sin(phi)
        row['shower_theta'] = event['shower_theta']
        row['shower_phi'] = event['shower_phi']
        row['shower_size'] = event['shower_size']
        row.append()


class KascadeLdfSimulation(BaseLdfSimulation):
    def __init__(self, *args, **kwargs):
        Ne = kwargs.pop('Ne', None)
        s = kwargs.pop('s', None)

        super(KascadeLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = KascadeLdf(Ne, s)

    def calculate_ldf_value(self, r, shower_size):
        return self.ldf.get_ldf_value_for_size(r, shower_size)


class BaseLdf(object):
    def calculate_core_distance_from_coordinates_and_direction(self,
                                                               x0, y0, x1, y1,
                                                               theta, phi):
        """Calculate core distance

        The core distance is the distance of the detector to the shower core,
        measured *on the shower front*.  For derivations, see logbook.

        """
        x = x1 - x0
        y = y1 - y0

        return sqrt(x ** 2 + y ** 2 -
                    (x * cos(phi) + y * sin(phi)) ** 2 * sin(theta) ** 2)


class KascadeLdf(BaseLdf):
    # shower parameters
    _Ne = 10 ** 4.8
    _s = .94
    _r0 = 40.
    _alpha = 1.5
    _beta = 3.6

    def __init__(self, Ne=None, s=None):
        if Ne is not None:
            self._Ne = Ne
        if s is not None:
            self._s = s

        self._cache_c_s_value()

    def _cache_c_s_value(self):
        self._c_s = self._c(self._s)

    def calculate_ldf_value(self, r):
        return self.get_ldf_value_for_size_and_shape(r, self._Ne, self._s)

    def get_ldf_value_for_size(self, r, Ne):
        return self.get_ldf_value_for_size_and_shape(r, Ne, self._s)

    def get_ldf_value_for_size_and_shape(self, r, Ne, s):
        c_s = self._c_s
        r0 = self._r0
        alpha = self._alpha
        beta = self._beta

        return Ne * c_s * (r / r0) ** (s - alpha) * (1 + r / r0) ** (s - beta)

    def _c(self, s):
        r0 = self._r0
        beta = self._beta
        alpha = self._alpha
        return gamma(beta - s) / (2 * pi * r0 ** 2 * gamma(s - alpha + 2) * gamma(alpha + beta - 2 * s - 2))
