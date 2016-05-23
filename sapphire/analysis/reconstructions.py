from itertools import izip_longest
import os

from numpy import isnan, histogram, linspace, percentile, std
from scipy.optimize import curve_fit
import tables

from ..storage import ReconstructedEvent, ReconstructedCoincidence
from ..clusters import HiSPARCStations, Station
from .. import api
from .direction_reconstruction import (EventDirectionReconstruction,
                                       CoincidenceDirectionReconstruction)
from .core_reconstruction import (EventCoreReconstruction,
                                  CoincidenceCoreReconstruction)
from .coincidence_queries import CoincidenceQuery
from .calibration import determine_detector_timing_offsets
from .event_utils import station_arrival_time
from ..utils import pbar, gauss, c


class ReconstructESDEvents(object):

    """Reconstruct events from single stations

    Example usage::

        >>> import tables
        >>> from sapphire import ReconstructESDEvents

        >>> data = tables.open_file('2014_1_1.h5', 'a')
        >>> station_path = '/hisparc/cluster_amsterdam/station_506'
        >>> rec = ReconstructESDEvents(data, station_path, 506, overwrite=True)
        >>> rec.reconstruct_and_store()

    To visualize the results::

        >>> import matplotlib.pyplot as plt
        >>> plt.polar([p for p in rec.phi if not isnan(p)],
        ...           [t for t in rec.theta if not isnan(t)], 'ko', alpha=0.2)

    or::

        >>> plt.polar(rec.reconstructions.col('azimuth'),
        ...           rec.reconstructions.col('zenith'), 'ko', alpha=0.2)

    """

    def __init__(self, data, station_group, station,
                 overwrite=False, progress=True,
                 destination='reconstructions'):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param station_group: the group containing the event table,
            the results will also be stored in this group.
        :param station: either a station number or
            :class:`~sapphire.clusters.Station` object. If number the
            positions and offsets are retrieved from the API. Otherwise
            the offsets will be determined with the available data.
        :param overwrite: if True, overwrite existing reconstruction table.
        :param progress: if True, show a progressbar while reconstructing.
        :param destination: alternative name for reconstruction table.

        """
        self.data = data
        self.station_group = data.get_node(station_group)
        self.events = self.station_group.events
        self.overwrite = overwrite
        self.progress = progress
        self.destination = destination
        self.offsets = [0., 0., 0., 0.]

        if isinstance(station, Station):
            self.station = station
            self.api_station = None
        else:
            cluster = HiSPARCStations([station])
            self.station = cluster.get_station(station)
            self.api_station = api.Station(station)

        self.direction = EventDirectionReconstruction(self.station)
        self.core = EventCoreReconstruction(self.station)

        self.theta = []
        self.phi = []
        self.detector_ids = []
        self.core_x = []
        self.core_y = []

    def reconstruct_and_store(self, detector_ids=None):
        """Shorthand function to reconstruct event and store the results"""

        self.prepare_output()
        if self.api_station is None:
            self.offsets = determine_detector_timing_offsets(self.events,
                                                             self.station)
            self.store_offsets()
        else:
            self.offsets = self.api_station
        self.reconstruct_directions(detector_ids=detector_ids)
        self.reconstruct_cores(detector_ids=detector_ids)
        self.store_reconstructions()

    def reconstruct_directions(self, detector_ids=None):
        """Reconstruct direction for all events

        :param detector_ids: list of detector ids to use for reconstructions.

        """
        angles = self.direction.reconstruct_events(self.events, detector_ids,
                                                   self.offsets, self.progress)
        self.theta, self.phi, self.detector_ids = angles

    def reconstruct_cores(self, detector_ids=None):
        """Reconstruct core for all events

        :param detector_ids: list of detector ids to use for reconstructions.

        """
        cores = self.core.reconstruct_events(self.events, detector_ids,
                                             self.progress)
        self.core_x, self.core_y = cores

    def prepare_output(self):
        """Prepare output table"""

        if self.destination in self.station_group:
            if self.overwrite:
                self.data.remove_node(self.station_group, self.destination,
                                      recursive=True)
            else:
                raise RuntimeError("Reconstructions table already exists for "
                                   "%s, and overwrite is False" %
                                   self.station_group)
        self.reconstructions = self.data.create_table(
            self.station_group, self.destination, ReconstructedEvent,
            expectedrows=self.events.nrows)
        self.reconstructions._v_attrs.station = self.station

    def store_offsets(self):
        """Store the determined offset in a table."""

        if 'detector_offsets' in self.station_group:
            if self.overwrite:
                self.data.remove_node(self.station_group.detector_offsets,
                                      recursive=True)
            else:
                raise RuntimeError("Detector offset table already exists for "
                                   "%s, and overwrite is False" %
                                   self.station_group)
        self.detector_offsets = self.data.create_array(
            self.station_group, 'detector_offsets', self.offsets)
        self.detector_offsets.flush()

    def store_reconstructions(self):
        """Loop over list of reconstructed data and store results

        Unsuccessful reconstructions are also stored but with the NumPy
        NaN as reconstructed value.

        """
        for event, core_x, core_y, theta, phi, detector_ids in izip_longest(
                self.events, self.core_x, self.core_y,
                self.theta, self.phi, self.detector_ids):
            self._store_reconstruction(event, core_x, core_y, theta, phi,
                                       detector_ids)
        self.reconstructions.flush()

    def _store_reconstruction(self, event, core_x, core_y, theta, phi,
                              detector_ids):
        """Store single reconstruction"""

        row = self.reconstructions.row
        row['id'] = event['event_id']
        row['ext_timestamp'] = event['ext_timestamp']
        try:
            row['min_n'] = min([event['n%d' % (id + 1)] for id in
                                detector_ids])
        except ValueError:
            # sometimes, all arrival times are -999 or -1, and then
            # detector_ids = []. So min([]) gives a ValueError.
            row['min_n'] = -999.
        row['x'] = core_x
        row['y'] = core_y
        row['zenith'] = theta
        row['azimuth'] = phi
        for id in detector_ids:
            row['d%d' % (id + 1)] = True
        row.append()


class ReconstructESDEventsFromSource(ReconstructESDEvents):

    def __init__(self, source_data, dest_data, source_group, dest_group,
                 station, overwrite=False, progress=True,
                 destination='reconstructions'):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param station_group: the group containing the event table,
            the results will also be stored in this group.
        :param station: either a station number or
            :class:`~sapphire.clusters.Station` object. If number the
            positions and offsets are retrieved from the API. Otherwise
            the offsets will be determined with the available data.
        :param overwrite: if True, overwrite existing reconstruction table.
        :param progress: if True, show a progressbar while reconstructing.
        :param destination: alternative name for reconstruction table.

        """
        super(ReconstructESDEventsFromSource, self).__init__(
            source_data, source_group, station, overwrite, progress,
            destination)
        self.dest_data = dest_data
        self.dest_group = dest_group

    def prepare_output(self):
        """Prepare output table"""

        dest_path = os.path.join(self.dest_group, self.destination)

        if dest_path in self.dest_data:
            if self.overwrite:
                self.dest_data.remove_node(dest_path, recursive=True)
            else:
                raise RuntimeError("Reconstructions table already exists for "
                                   "%s, and overwrite is False" %
                                   self.dest_group)
        self.reconstructions = self.dest_data.create_table(
            self.dest_group, self.destination, ReconstructedEvent,
            expectedrows=self.events.nrows, createparents=True)
        self.reconstructions._v_attrs.station = self.station


class ReconstructESDCoincidences(object):

    """Reconstruct coincidences, e.g. event between multiple stations

    Example usage::

        >>> import tables
        >>> from sapphire import ReconstructESDCoincidences

        >>> data = tables.open_file('2014_1_1.h5', 'a')
        >>> rec = ReconstructESDCoincidences(data, overwrite=True)
        >>> rec.reconstruct_and_store()

    """

    def __init__(self, data, coincidences_group='/coincidences',
                 overwrite=False, progress=True,
                 destination='reconstructions', cluster=None):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param coincidences_group: the destination group.
        :param overwrite: if True, overwrite existing reconstruction table.
        :param progress: if True, show a progressbar while reconstructing.
        :param destination: alternative name for reconstruction table.
        :param cluster: a Cluster object to use for the reconstructions.

        """
        self.data = data
        self.coincidences_group = data.get_node(coincidences_group)
        self.coincidences = self.coincidences_group.coincidences
        self.overwrite = overwrite
        self.progress = progress
        self.destination = destination
        self.offsets = {}

        self.cq = CoincidenceQuery(data, self.coincidences_group)
        if cluster is None:
            try:
                self.cluster = self.coincidences_group._f_getattr('cluster')
            except AttributeError:
                s_active = self._get_active_stations()
                self.cluster = HiSPARCStations(s_active)
        else:
            self.cluster = cluster

        self.direction = CoincidenceDirectionReconstruction(self.cluster)
        self.core = CoincidenceCoreReconstruction(self.cluster)

        self.theta = []
        self.phi = []
        self.station_numbers = []
        self.core_x = []
        self.core_y = []

    def reconstruct_and_store(self, station_numbers=None):
        """Shorthand function to reconstruct coincidences and store results"""

        self.prepare_output()
        self.get_station_timing_offsets()
        self.reconstruct_directions(station_numbers=station_numbers)
        self.reconstruct_cores(station_numbers=station_numbers)
        self.store_reconstructions()

    def reconstruct_directions(self, station_numbers=None):
        """Reconstruct direction for all events

        :param detector_ids: list of detector ids to use for reconstructions.

        """
        coincidences = pbar(self.cq.all_coincidences(iterator=True),
                            length=self.coincidences.nrows, show=self.progress)
        angles = self.direction.reconstruct_coincidences(
            self.cq.all_events(coincidences, n=0), station_numbers,
            self.offsets, progress=False)
        self.theta, self.phi, self.station_numbers = angles

    def reconstruct_cores(self, station_numbers=None):
        """Reconstruct core for all events

        :param detector_ids: list of detector ids to use for reconstructions.

        """
        coincidences = pbar(self.cq.all_coincidences(iterator=True),
                            length=self.coincidences.nrows, show=self.progress)
        cores = self.core.reconstruct_coincidences(
            self.cq.all_events(coincidences, n=0), station_numbers,
            progress=False)
        self.core_x, self.core_y = cores

    def prepare_output(self):
        """Prepare output table"""

        if self.destination in self.coincidences_group:
            if self.overwrite:
                self.data.remove_node(self.coincidences_group,
                                      self.destination, recursive=True)
            else:
                raise RuntimeError("Reconstructions table already exists for "
                                   "%s, and overwrite is False" %
                                   self.coincidences_group)

        s_columns = {'s%d' % station.number: tables.BoolCol(pos=p)
                     for p, station in enumerate(self.cluster.stations, 26)}
        description = ReconstructedCoincidence
        description.columns.update(s_columns)
        self.reconstructions = self.data.create_table(
            self.coincidences_group, self.destination, description,
            expectedrows=self.coincidences.nrows)
        self.reconstructions._v_attrs.cluster = self.cluster

    def get_station_timing_offsets(self):
        """Get predetermined station offsets

        References for subclusters:
        102 for Zaanlands stations, data from 2012/6-2014/8.
        501 for Science Park stations, data from 2010/1-2014/8.
        7001 for Twente University stations, data from 2011/8-2014/8.
        8001 for Eindhoven University stations, data from 2011/10-2014/8.

        """
        self.offsets = {102: [-3.1832, 0.0000, 0.0000, 0.0000],
                        104: [-1.5925, -5.0107, 0.0000, 0.0000],
                        105: [-14.1325, -10.9451, 0.0000, 0.0000],
                        501: [-1.10338, 0.0000, 5.35711, 3.1686],
                        502: [-8.11711, -8.5528, -8.72451, -9.3388],
                        503: [-22.9796, -26.6098, -22.7522, -21.8723],
                        504: [-15.4349, -15.2281, -15.1860, -16.5545],
                        505: [-21.6035, -21.3060, -19.6826, -25.5366],
                        506: [-20.2320, -15.8309, -14.1818, -14.1548],
                        508: [-26.2402, -24.9859, -24.0131, -23.2882],
                        509: [-24.8369, -23.0218, -20.6011, -24.3757],
                        7001: [4.5735, 0.0000, 0.0000, 0.0000],
                        7002: [45.0696, 47.8311, 0.0000, 0.0000],
                        7003: [-2.2674, -4.9578, 0.0000, 0.0000],
                        8001: [2.5733, 0.0000, 0.0000, 0.0000],
                        8004: [-39.3838, -36.1131, 0.0000, 0.0000],
                        8008: [57.3990, 58.1135, 0.0000, 0.0000],
                        8009: [-20.3489, -16.9938, 0.0000, 0.0000]}

    def determine_station_timing_offsets(self, ref_station_number=501):
        """Determine the offsets between the stations.

        ADL: This should use more than one day of data for a good fit.
        Station altitudes are not taken into account. It would be better
        to choose a different reference station per (sub)cluster.

        """
        # First determine detector offsets for each station
        for s_path in self.coincidences_group.s_index:
            try:
                station_group = self.data.get_node(s_path)
            except tables.NoSuchNodeError:
                continue
            station_number = int(s_path.split('station_')[-1])
            station = self.cluster.get_station(station_number)
            offsets = determine_detector_timing_offsets(station_group.events,
                                                        station)
            self.offsets[station_number] = offsets

        # Currently disabled station offsets because they do not work well.

        # Now determine station offsets and add those to detector offsets
        ref_station = self.cluster.get_station(ref_station_number)
        ref_d_ids = range(len(ref_station.detectors))
        ref_z = ref_station.calc_center_of_mass_coordinates()[2]
        ref_d_off = self.offsets[ref_station_number]

        for station in pbar(self.cluster.stations):
            # Skip reference station
            if station.number == ref_station_number:
                continue
            d_ids = range(len(station.detectors))
            z = station.calc_center_of_mass_coordinates()[2]
            d_off = self.offsets[station.number]

            stations = [ref_station_number, station.number]
            coincidences = self.cq.all(stations)
            c_events = self.cq.events_from_stations(coincidences, stations)

            dt = []
            for events in c_events:
                # Filter for possibility of same station twice in coincidence
                if len(events) is not 2:
                    continue
                if events[0][0] == ref_station_number:
                    ref_event = events[0][1]
                    event = events[1][1]
                else:
                    ref_event = events[1][1]
                    event = events[0][1]

                ref_ts = ref_event['ext_timestamp']
                ref_t = station_arrival_time(ref_event, ref_ts, ref_d_ids,
                                             ref_d_off)
                if isnan(ref_t):
                    continue
                t = station_arrival_time(event, ref_ts, d_ids, d_off)
                if isnan(t):
                    continue
                dt.append(t - ref_t)

            bins = linspace(percentile(dt, 2), percentile(dt, 98), 200)
            y, bins = histogram(dt, bins=bins)
            x = (bins[:-1] + bins[1:]) / 2
            try:
                popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0., std(dt)))
                station_offset = popt[1] + (z - ref_z) / c
            except RuntimeError:
                station_offset = 0.
            self.offsets[station.number] = [detector_offset + station_offset
                                            for detector_offset in
                                            self.offsets[station.number]]

    def store_reconstructions(self):
        """Loop over list of reconstructed data and store results

        Unsuccessful reconstructions are also stored but with the NumPy
        NaN as reconstructed value.

        """
        for coincidence, x, y, theta, phi, station_numbers in izip_longest(
                self.coincidences, self.core_x, self.core_y,
                self.theta, self.phi, self.station_numbers):
            self._store_reconstruction(coincidence, x, y, theta, phi,
                                       station_numbers)
        self.reconstructions.flush()

    def _store_reconstruction(self, coincidence, core_x, core_y, theta, phi,
                              station_numbers):
        """Store single reconstruction"""

        row = self.reconstructions.row

        row['id'] = coincidence['id']
        row['ext_timestamp'] = coincidence['ext_timestamp']
        row['x'] = core_x
        row['y'] = core_y
        row['zenith'] = theta
        row['azimuth'] = phi

        row['reference_x'] = coincidence['x']
        row['reference_y'] = coincidence['y']
        row['reference_zenith'] = coincidence['zenith']
        row['reference_azimuth'] = coincidence['azimuth']
        row['reference_size'] = coincidence['size']
        row['reference_energy'] = coincidence['energy']

        for number in station_numbers:
            row['s%d' % number] = True

        row.append()

    def _get_active_stations(self):
        """Return station numbers with non-empty event table in datafile"""

        active_stations = []

        for s_path in self.coincidences_group.s_index:
            try:
                station_event_table = self.data.get_node(s_path + '/events')
            except tables.NoSuchNodeError:
                continue
            if not station_event_table.nrows:
                continue
            active_stations.append(int(s_path.split('station_')[-1]))

        return active_stations


class ReconstructESDCoincidencesFromSource(ReconstructESDCoincidences):

    def __init__(self, source_data, dest_data, source_group, dest_group,
                 overwrite=False, progress=True,
                 destination='reconstructions', cluster=None):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param station_group: the group containing the event table,
            the results will also be stored in this group.
        :param station: either a station number or
            :class:`~sapphire.clusters.Station` object. If number the
            positions and offsets are retrieved from the API. Otherwise
            the offsets will be determined with the available data.
        :param overwrite: if True, overwrite existing reconstruction table.
        :param progress: if True, show a progressbar while reconstructing.
        :param destination: alternative name for reconstruction table.

        """
        super(ReconstructESDCoincidencesFromSource, self).__init__(
            source_data, source_group, overwrite, progress, destination,
            cluster)
        self.dest_data = dest_data
        self.dest_group = dest_group

    def prepare_output(self):
        """Prepare output table"""

        dest_path = os.path.join(self.dest_group, self.destination)

        if dest_path in self.dest_data:
            if self.overwrite:
                self.dest_data.remove_node(dest_path, recursive=True)
            else:
                raise RuntimeError("Reconstructions table already exists for "
                                   "%s, and overwrite is False" %
                                   self.dest_group)

        s_columns = {'s%d' % station.number: tables.BoolCol(pos=p)
                     for p, station in enumerate(self.cluster.stations, 26)}
        description = ReconstructedCoincidence
        description.columns.update(s_columns)
        self.reconstructions = self.dest_data.create_table(
            self.dest_group, self.destination, description,
            expectedrows=self.coincidences.nrows, createparents=True)
        self.reconstructions._v_attrs.cluster = self.cluster
