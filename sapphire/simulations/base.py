class BaseSimulation(object):

    """Base class for simulations.

    :param cluster: :class:`sapphire.clusters.BaseCluster` instance.
    :param datafile: writeable PyTables file handle.
    :param output_path: path (as string) to the PyTables group (need not
                        exist) in which the result tables will be created.
    :param N: number of simulations to perform.
    
    """

    def __init__(self, cluster, datafile, output_path='/', N=1):
        self.cluster = cluster
        self.datafile = datafile
        self.output_path = output_path
        self.N = N

        self._prepare_output_tables()

    def _prepare_output_tables(self):
        """Prepare output tables in datafile."""

        pass

    def run(self):
        """Run the simulations."""

        for shower_id in range(self.N):
            station_events = []
            shower_parameters = self.generate_shower_parameters()

            for station in self.cluster.stations:
                has_triggered, station_observables = \
                    self.simulate_station_response(station,
                                                   shower_parameters)
                if has_triggered:
                    self.store_station_observables(station,
                                                   station_observables)
                    station_events.append((station, station_observables))

            self.store_coincidence(shower_id, shower_parameters,
                                   station_events)

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc."""

        pass

    def simulate_station_response(self, station, shower_parameters):
        """Simulate station response to a shower."""

        station_observables = []
        for detector in station.detectors:
            observables = self.simulate_detector_response(detector,
                                                          shower_parameters)
            station_observables.append(observables)

        has_triggered = self.simulate_trigger(station_observables)
        return has_triggered, station_observables
    
    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower."""

        # implement this!
        observables = None

        return observables

    def simulate_trigger(self, station_observables):
        """Simulate a trigger response."""

        return True

    def store_station_observables(self, station, station_observables):
        """Store station observables."""

        pass

    def store_coincidence(self, shower_id, shower_parameters,
                          station_events):
        """Store coincidence."""

        pass
