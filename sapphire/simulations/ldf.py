from math import pi, sin, cos, sqrt
from scipy.special import gamma
import progressbar as pb
import sys

from base import BaseSimulation


class BaseLdfSimulation(BaseSimulation):
    """ Simulation using a lateral distribution function as a model

        Not using an EAS simulation but rather a lateral distribution function
        as the model for the detector simulation considerably speeds up the
        complete process.  However, we will have to validate the use of an LDF.

        This is the base class, so there is no LDF defined.  Override
        calculate_ldf_value to complete an implementation.

    """
    # Number of MIPs required to go over threshold
    trig_threshold = .8

    # Alternative for observables.nrows, which is only updated with a flush()
    # We don't want to flush often, so we use (and update!) this instead
    _observables_nrows = 0


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
                     'shower_theta': 0., 'shower_phi': 0.}
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

        num_detectors_over_threshold = sum([True if u >= 1 else False for u
                                            in num_particles])
        if num_detectors_over_threshold >= 2:
            has_triggered = True
        else:
            has_triggered = False

        return has_triggered, id

    def simulate_detector_observables(self, detector, event):
        R = self.calculate_core_distance(detector, event)
        N = self.calculate_ldf_value(R)
        return N * detector.get_area()

    def calculate_core_distance(self, detector, event):
        r, phi = event['r'], event['phi']
        x = r * cos(phi)
        y = r * sin(phi)

        X, Y = detector.get_xy_coordinates()

        return sqrt((x - X) ** 2 + (y - Y) ** 2)

    def calculate_ldf_value(self, R):
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


class KascadeLdfSimulation(BaseLdfSimulation):
    # shower parameters
    _Ne = 10 ** 4.8
    _s = .94
    _r0 = 40.
    _alpha = 1.5
    _beta = 3.6

    def __init__(self, *args, **kwargs):
        self.cache_c_s_value()

        super(KascadeLdfSimulation, self).__init__(*args, **kwargs)

    def cache_c_s_value(self):
        self._c_s = self._c(self._s)

    def run(self):
        self.cache_c_s_value()

        super(KascadeLdfSimulation, self).run()

    def calculate_ldf_value(self, r):
        Ne = self._Ne
        s = self._s
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
