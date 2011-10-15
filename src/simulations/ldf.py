from math import sin, cos, sqrt

from base import BaseSimulation

class BaseLdfSimulation(BaseSimulation):
    """ Simulation using a lateral distribution function as a model

        Not using an EAS simulation but rather a lateral distribution function
        as the model for the detector simulation considerably speeds up the
        complete process.  However, we will have to validate the use of an LDF.

        This is the base class, so there is no LDF defined.  Override
        calculate_ldf_value to complete an implementation.

    """
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

        for event_id, (r, phi) in enumerate(positions):
            event = {'id': event_id, 'r': r, 'phi': phi}
            self.simulate_event(event)

        self._run_exit_msg()

    def simulate_event(self, event):
        multiplicity = 0
        station_event_ids = []

        for station in self.cluster.stations:
            has_trigger, id = self.simulate_station_observables_and_return_has_triggered_and_eventid(station, event)
            multiplicity += has_trigger
            station_event_ids.append(id)
        self.write_coincidence(event, multiplicity)
        self.c_index.append((event['id'], station_event_ids))

    def simulate_station_observables_and_return_has_triggered_and_eventid(self, station, event):
        num_particles = []
        for detector in station.detectors:
            num_particles.append(self.simulate_detector_observables(detector, event))
        self.write_observables(station, *num_particles)

        return sum([True if u >= 1 else False for u in num_particles])

    def simulate_detector_observables(self, detector, event):
        R = self.calculate_core_distance(detector, event)
        N = self.calculate_ldf_value(R)
        return N * detector.get_area()

    def calculate_core_distance(self, detector, event):
        r, phi = event['r'], event['phi']
        x = r * cos(phi)
        y = r * sin(phi)

        X, Y = detector.get_position()

        return sqrt((x - X) ** 2 + (y - Y) ** 2)

    def calculate_ldf_value(self, R):
        return 0.

    def write_observables(self, station, event, n1, n2, n3, n4):
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
        row['N'] = sum([n1, n2, n3, n4])
        row['t1'], row['t2'], row['t3'], row['t4'] = 0, 0, 0, 0
        row['n1'], row['n2'], row['n3'], row['n4'] = n1, n2, n3, n4
        row.append()
