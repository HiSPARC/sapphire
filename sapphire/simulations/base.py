"""Perform simulations of air showers on a cluster of stations

This base class can be subclassed to provide various kinds of
simulations. These simulations will inherit the base functionallity from
this class, including the creation of event and coincidence tables to
store the results, which will look similar to regular HiSPARC data, such
that the same reconstruction analysis can be applied to both.

Example usage::

    >>> import tables

    >>> from sapphire.simulations.base import BaseSimulation
    >>> from sapphire import ScienceParkCluster

    >>> data = tables.open_file('/tmp/test_base_simulation.h5', 'w')
    >>> cluster = ScienceParkCluster()

    >>> sim = BaseSimulation(cluster, data, '/simulations/this_run', 10)
    >>> sim.run()

"""
import warnings
import random

import numpy as np
import tables

from .. import storage
from ..analysis.process_events import ProcessEvents
from ..utils import pbar


class BaseSimulation(object):

    """Base class for simulations.

    :param cluster: :class:`~sapphire.clusters.BaseCluster` instance.
    :param data: writeable PyTables file handle.
    :param output_path: path (as string) to the PyTables group (need not
                        exist) in which the result tables will be created.
    :param N: number of simulations to perform.
    :param seed: seed for the pseudo-random number generators.
    :param progress: if True, show a progressbar while simulating.

    """

    def __init__(self, cluster, data, output_path='/', N=1, seed=None,
                 progress=True):
        self.cluster = cluster
        self.data = data
        self.output_path = output_path
        self.N = N
        self.progress = progress

        self._prepare_output_tables()

        if seed is not None:
            random.seed(seed)
            np.random.seed(seed)

    def _prepare_output_tables(self):
        """Prepare output tables in the output data file.

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

        for (shower_id, shower_parameters) in enumerate(
                self.generate_shower_parameters()):

            station_events = self.simulate_events_for_shower(shower_parameters)
            self.store_coincidence(shower_id, shower_parameters,
                                   station_events)

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc."""

        shower_parameters = {'core_pos': (None, None),
                             'zenith': None,
                             'azimuth': None,
                             'size': None,
                             'energy': None,
                             'ext_timestamp': None}

        for i in pbar(range(self.N), show=self.progress):
            yield shower_parameters

    def simulate_events_for_shower(self, shower_parameters):
        """Simulate station events for a single shower"""

        station_events = []
        for station_id, station in enumerate(self.cluster.stations):
            has_triggered, station_observables = \
                self.simulate_station_response(station,
                                               shower_parameters)
            if has_triggered:
                event_index = \
                    self.store_station_observables(station_id,
                                                   station_observables)
                station_events.append((station_id, event_index))
        return station_events

    def simulate_station_response(self, station, shower_parameters):
        """Simulate station response to a shower."""

        detector_observables = self.simulate_all_detectors(
            station.detectors, shower_parameters)
        has_triggered = self.simulate_trigger(detector_observables)
        station_observables = \
            self.process_detector_observables(detector_observables)
        station_observables = self.simulate_gps(station_observables,
                                                shower_parameters, station)

        return has_triggered, station_observables

    def simulate_all_detectors(self, detectors, shower_parameters):
        """Simulate response of all detectors in a station.

        :param detectors: list of detectors
        :param shower_parameters: parameters of the shower

        """
        detector_observables = []
        for detector in detectors:
            observables = self.simulate_detector_response(detector,
                                                          shower_parameters)
            detector_observables.append(observables)

        return detector_observables

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        :param detector: :class:`~sapphire.clusters.Detector` instance
        :param shower_parameters: shower parameters
        :return: dictionary with keys 'n' (number of particles in
            detector) and 't' (time of arrival of first detected particle).

        """
        # implement this!
        observables = {'n': 0., 't': -999}

        return observables

    def simulate_trigger(self, detector_observables):
        """Simulate a trigger response."""

        return True

    def simulate_gps(self, station_observables, shower_parameters, station):
        """Simulate gps timestamp."""

        gps_timestamp = {'ext_timestamp': 0, 'timestamp': 0, 'nanoseconds': 0}
        station_observables.update(gps_timestamp)

        return station_observables

    def process_detector_observables(self, detector_observables):
        """Process detector observables for a station.

        The list of detector observables is converted into a dictionary
        containing the familiar observables like pulseheights, n1, n2,
        ..., t1, t2, ..., integrals, etc.

        :param detector_observables: list of observables of the detectors
                                     making up a station.
        :return: dictionary containing the familiar station observables
                 like n1, n2, n3, etc.

        """
        station_observables = {'pulseheights': 4 * [-1.],
                               'integrals': 4 * [-1.]}

        for detector_id, observables in enumerate(detector_observables, 1):
            for key, value in observables.iteritems():
                if key in ['n', 't']:
                    key = key + str(detector_id)
                    station_observables[key] = value
                elif key in ['pulseheights', 'integrals']:
                    idx = detector_id - 1
                    station_observables[key][idx] = value

        return station_observables

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
        for key, value in station_observables.iteritems():
            if key in events_table.colnames:
                row[key] = value
            else:
                warnings.warn('Unsupported variable')
        row.append()
        events_table.flush()

        return events_table.nrows - 1

    def store_coincidence(self, shower_id, shower_parameters,
                          station_events):
        """Store coincidence.

        Store the information to find events of different stations
        belonging to the same simulated shower in the coincidences
        tables.

        :param shower_id: The shower number for the coincidence id.
        :param shower_parameters: A dictionary with the parameters of
            the simulated shower.
        :param station_events: A list of tuples containing the
            station_id and event_index referring to the events that
            participated in the coincidence.

        """
        row = self.coincidences.row
        row['id'] = shower_id
        row['N'] = len(station_events)
        row['x'], row['y'] = shower_parameters['core_pos']
        row['zenith'] = shower_parameters['zenith']
        row['azimuth'] = shower_parameters['azimuth']
        row['size'] = shower_parameters['size']
        row['energy'] = shower_parameters['energy']

        timestamps = []
        for station_id, event_index in station_events:
            station = self.cluster.stations[station_id]
            row['s%d' % station.number] = True
            station_group = self.station_groups[station_id]
            event = station_group.events[event_index]
            timestamps.append((event['ext_timestamp'], event['timestamp'],
                               event['nanoseconds']))

        try:
            first_timestamp = sorted(timestamps)[0]
        except IndexError:
            first_timestamp = (0, 0, 0)

        row['ext_timestamp'], row['timestamp'], row['nanoseconds'] = \
            first_timestamp
        row.append()
        self.coincidences.flush()

        self.c_index.append(station_events)
        self.c_index.flush()

    def _prepare_coincidence_tables(self):
        """Create coincidence tables

        These are the same as the tables created by
        :class:`~sapphire.analysis.coincidences.CoincidencesESD`.
        This makes it easy to link events detected by multiple stations.

        """
        self.coincidence_group = self.data.create_group(self.output_path,
                                                        'coincidences',
                                                        createparents=True)
        self.coincidence_group._v_attrs.cluster = self.cluster

        description = storage.Coincidence
        s_columns = {'s%d' % station.number: tables.BoolCol(pos=p)
                     for p, station in enumerate(self.cluster.stations, 12)}
        description.columns.update(s_columns)

        self.coincidences = self.data.create_table(
            self.coincidence_group, 'coincidences', description)

        self.c_index = self.data.create_vlarray(
            self.coincidence_group, 'c_index', tables.UInt32Col(shape=2))

        self.s_index = self.data.create_vlarray(
            self.coincidence_group, 's_index', tables.VLStringAtom())

    def _prepare_station_tables(self):
        """Create the groups and events table to store the observables

        :param id: the station number, used for the group name
        :param station: a :class:`~sapphire.clusters.Station` object

        """
        self.cluster_group = self.data.create_group(self.output_path,
                                                    'cluster_simulations',
                                                    createparents=True)
        self.station_groups = []
        for station in self.cluster.stations:
            station_group = self.data.create_group(self.cluster_group,
                                                   'station_%d' %
                                                   station.number)
            description = ProcessEvents.processed_events_description
            self.data.create_table(station_group, 'events', description,
                                   expectedrows=self.N)
            self.station_groups.append(station_group)

    def _store_station_index(self):
        """Stores the references to the station groups for coincidences"""

        for station_group in self.station_groups:
            self.s_index.append(station_group._v_pathname)
        self.s_index.flush()
