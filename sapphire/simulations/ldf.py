from math import pi, sin, cos, sqrt
from scipy.special import gamma
import progressbar as pb
import sys

import numpy as np

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


    def __init__(self, cluster, data, output, R, N, shower_size=10 ** 4.8, use_poisson=None, gauss=None, trig_threshold=1., **kwargs):
        self.use_poisson = use_poisson
        self.gauss = gauss
        self.trig_threshold = trig_threshold
        self.shower_size = shower_size

        super(BaseLdfSimulation, self).__init__(cluster, data, output, R, N, **kwargs)

    def run(self, positions=None):
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
            event = {'id': event_id, 'r': r, 'phi': phi, 'alpha': 0.,
                     'shower_theta': 0., 'shower_phi': 0.,
                     'shower_size': self.shower_size}
            self.simulate_event(event)

        self.coincidences.flush()
        self.observables.flush()
        # integrity check
        assert(self._observables_nrows == self.observables.nrows)
        self.c_index.flush()

        self._run_exit_msg()

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
        id = self.write_observables_and_return_id(station, event, *num_particles)

        num_detectors_over_threshold = sum([True if u >= self.trig_threshold else False for u
                                            in num_particles])
        if num_detectors_over_threshold >= 2:
            has_triggered = True
        else:
            has_triggered = False

        return has_triggered, id

    def simulate_detector_observables(self, detector, event):
        R = self.calculate_core_distance(detector, event)
        density = self.calculate_ldf_value(R, event['shower_size'])
        num_particles = density * detector.get_area()

        if self.use_poisson:
            N = np.random.poisson(num_particles)
        else:
            N = num_particles

        if self.gauss is not None and N > 0:
            N = np.random.normal(loc=N, scale=sqrt(N) * self.gauss)

        return N

    def calculate_core_distance(self, detector, event):
        r, phi = event['r'], event['phi']
        x = r * cos(phi)
        y = r * sin(phi)

        X, Y = detector.get_xy_coordinates()

        return sqrt((x - X) ** 2 + (y - Y) ** 2)

    def calculate_ldf_value(self, R, shower_size):
        return 0.

    def write_observables_and_return_id(self, station, event, n1, n2, n3, n4):
        """Write observables from a single event

        :param station: Station instance
        :param n1, n2, n3, n4: number of particles detected for each detector
            1, 2, 3 and 4

        """
        row = self.observables.row

        x, y, alpha = station.get_xyalpha_coordinates()
        r, phi, alpha = station.get_rphialpha_coordinates()

        row['id'] = event['id']
        row['station_id'] = station.station_id
        row['r'] = r
        row['phi'] = phi
        row['x'] = x
        row['y'] = y
        row['alpha'] = alpha
        row['N'] = sum([1 if u >= self.trig_threshold else 0 for u
                        in n1, n2, n3, n4])
        row['t1'], row['t2'], row['t3'], row['t4'] = 0, 0, 0, 0
        row['n1'], row['n2'], row['n3'], row['n4'] = n1, n2, n3, n4
        row.append()
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


class KascadeLdf():
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
