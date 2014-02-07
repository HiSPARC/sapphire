import tables

from sapphire import storage


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
        """Prepare output tables in datafile.

        The groups and tables will be created in the output_path.

        :raises tables.NodeError: If any of the groups (e.g.
            '/coincidences') already exist a exception will be raised.
        :raises tables.FileModeError: If the datafile is not writeable.

        """
        self._prepare_coincidence_tables()
        self._prepare_station_tables()
        self._store_station_index()

    def run(self):
        """Run the simulations."""

        for shower_id in range(self.N):
            station_events = []
            shower_parameters = self.generate_shower_parameters()

            for station_id, station in enumerate(self.cluster.stations):
                has_triggered, station_observables = \
                        self.simulate_station_response(station,
                                                       shower_parameters)
                if has_triggered:
                    event_index = \
                            self.store_station_observables(station_id,
                                                           station_observables)
                    station_events.append((station_id, event_index))

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

    def store_station_observables(self, station_id, station_observables):
        """Store station observables.

        :param station_id: the id of the station in self.cluster
        :param station_observables: A dictionary containing the
            variables to be stored for this event.
        :return: The index (row number) of the newly added event.

        """
        events_table = self.station_groups[station_id].events
        row = events_table.row
        row['event_id'] = events_table.nrows
        # for key, value in station_observables.iteritems():
        #     if key in events_table.colnames:
        #         row[key] = value
        #     else:
        #         pass
        row.append()
        events_table.flush()

        return events_table.nrows - 1

    def store_coincidence(self, shower_id, shower_parameters,
                          station_events):
        """Store coincidence."""

        pass

    def _prepare_coincidence_tables(self):
        """Create coincidence tables

        These are the same as the tables created by
        :class:`sapphire.analysis.coincidences.CoincidencesESD`.
        This makes it easy to link events detected by multiple stations.

        """
        self.coincidence_group = self.datafile.createGroup(self.output_path,
                                                           'coincidences',
                                                           createparents=True)
        self.coincidence_group._v_attrs.cluster = self.cluster

        description = storage.Coincidence
        stations_description = {'s%d' % n: tables.BoolCol()
                                for n in range(len(self.cluster.stations))}
        description.columns.update(stations_description)

        self.coincidences = self.datafile.createTable(
                self.coincidence_group, 'coincidences', description)

        self.c_index = self.datafile.createVLArray(
                self.coincidence_group, 'c_index', tables.UInt32Col(shape=2))

        self.s_index = self.datafile.createVLArray(
                self.coincidence_group, 's_index', tables.VLStringAtom())

    def _prepare_station_tables(self):
        """Create the groups and events table to store the observables

        :param id: the station number, used for the group name
        :param station: a :class:`sapphire.clusters.Station` object

        """
        self.cluster_group = self.datafile.createGroup(self.output_path,
                                                      'cluster_simulations',
                                                      createparents=True)
        self.station_groups = []
        for id, station in enumerate(self.cluster.stations):
            station_group = self.datafile.createGroup(self.cluster_group,
                                                      'station_%d' % id)
            events_table = \
                    self.datafile.createTable(station_group, 'events',
                                              storage.ProcessedHisparcEvent,
                                              expectedrows=self.N)
            self.station_groups.append(station_group)

    def _store_station_index(self):
        """Stores the references to the station groups for coincidences"""

        for station_group in self.station_groups:
            self.s_index.append(station_group._v_pathname)
        self.s_index.flush()
