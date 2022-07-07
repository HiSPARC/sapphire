""" Reconstruct HiSPARC events and coincidences

    This module contains classes that can be used to reconstruct
    HiSPARC events and coincidences. These classes can be used to automate
    the tasks of reconstructing directions and/or cores.

    The classes can reconstruct measured data from the ESD as well as
    simulated data from :mod:`sapphire.simulations`.

    The classes read data stored in HDF5 files and extract station metadata
    (cluster and detector layout, station and detector offsets) from
    various sources:

    - from the public database using :class:`sapphire.api.Station` objects
    - from stored or provided :class`sappire.cluster.Station` objects,
      usually cluster or station layout stored by :mod:`sapphire.simulations`

     Reconstructed data is stored in HDF5 files.

"""
import os
import warnings

from itertools import zip_longest

import tables

from .. import api
from ..clusters import HiSPARCStations, Station
from ..storage import ReconstructedCoincidence, ReconstructedEvent
from ..utils import pbar
from .calibration import determine_detector_timing_offsets
from .coincidence_queries import CoincidenceQuery
from .core_reconstruction import CoincidenceCoreReconstruction, EventCoreReconstruction
from .direction_reconstruction import CoincidenceDirectionReconstruction, EventDirectionReconstruction


class ReconstructESDEvents:

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
                 overwrite=False, progress=True, verbose=False,
                 destination='reconstructions',
                 force_fresh=False, force_stale=False):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param station_group: the group containing the event table,
            the results will also be stored in this group.
        :param station: either a station number or
            :class:`sapphire.clusters.Station` object. If it is a number the
            positions and offsets will be retrieved from the public database
            or retrieved from the datafile when stored by a simulation.
            Otherwise the offsets will be determined with the available data.
        :param overwrite: if True overwrite existing reconstruction table.
        :param progress: if True show a progressbar while reconstructing.
        :param verbose: if True be verbose about station metadata usage.
        :param destination: alternative name for reconstruction table.

        """
        self.data = data
        self.station_group = data.get_node(station_group)
        self.events = self.station_group.events
        self.overwrite = overwrite
        self.progress = progress
        self.verbose = verbose
        self.destination = destination
        self.force_fresh = force_fresh
        self.force_stale = force_stale

        self.offsets = [0., 0., 0., 0.]

        self._get_or_create_station_object(station)

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
        self.get_detector_offsets()
        self.reconstruct_directions(detector_ids=detector_ids)
        self.reconstruct_cores(detector_ids=detector_ids)
        self.store_reconstructions()

    def reconstruct_directions(self, detector_ids=None):
        """Reconstruct direction for all events

        :param detector_ids: list of detector ids to use for reconstructions.

        """
        if len(self.core_x) and len(self.core_y):
            initials = ({'core_x': x, 'core_y': y}
                        for x, y in zip(self.core_x, self.core_y))
        else:
            initials = []
        angles = self.direction.reconstruct_events(self.events, detector_ids,
                                                   self.offsets, self.progress,
                                                   initials)
        self.theta, self.phi, self.detector_ids = angles

    def reconstruct_cores(self, detector_ids=None):
        """Reconstruct core for all events

        :param detector_ids: list of detector ids to use for reconstructions.

        """
        if len(self.theta) and len(self.phi):
            initials = ({'theta': theta, 'phi': phi}
                        for theta, phi in zip(self.theta, self.phi))
        else:
            initials = []
        cores = self.core.reconstruct_events(self.events, detector_ids,
                                             self.progress, initials)
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
        try:
            self.reconstructions._v_attrs.station = self.station
        except tables.HDF5ExtError:
            warnings.warn('Unable to store station object, to large for HDF.')

    def get_detector_offsets(self):
        """Get or determine detector offsets

        Try to extract the offsets from the provided cluster object.
        If those are not available use the :class:`sapphire.api.Station`
        object for the station number.
        Else determine the offsets from the event table.

        - if a cluster object is provided:

          -  use offsets from that object if available in the object
          -  else determine the offsets from the events in datafile with the
             provided cluster object.

        - if a station number is provided:

          -  if a cluster object is stored in the datafile use offsets from
             that object if available.
          -  else get offsets from `api.Station` object.

        """
        try:
            self.offsets = [d.offset for d in self.station.detectors]
            if self.verbose:
                print('Read detector offsets from station object.')
        except AttributeError:
            if self.station_number is not None:
                self.offsets = api.Station(self.station_number,
                                           force_fresh=self.force_fresh,
                                           force_stale=self.force_stale)
                if self.verbose:
                    print('Reading detector offsets from public database.')
            else:
                self.offsets = determine_detector_timing_offsets(self.events,
                                                                 self.station)
                self.store_offsets()
                if self.verbose:
                    print('Determined offsets from event data: ', self.offsets)

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
        for event, core_x, core_y, theta, phi, detector_ids in zip_longest(
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
            row['min_n'] = min(event['n%d' % (id + 1)] for id in detector_ids)
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

    def _get_or_create_station_object(self, station):
        if isinstance(station, Station):
            self.station = station
            self.station_number = None
            if self.verbose:
                print('Using object %s for metadata.' % self.station)
        else:
            self.station_number = station
            cluster = HiSPARCStations([station],
                                      force_fresh=self.force_fresh,
                                      force_stale=self.force_stale)
            self.station = cluster.get_station(station)
            if self.verbose:
                print(f'Constructed object {self.station} from public database.')


class ReconstructESDEventsFromSource(ReconstructESDEvents):

    def __init__(self, source_data, dest_data, source_group, dest_group,
                 station, overwrite=False, progress=True, verbose=False,
                 destination='reconstructions',
                 force_fresh=False, force_stale=False):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param station_group: the group containing the event table,
            the results will also be stored in this group.
        :param station: either a station number or
            :class:`sapphire.clusters.Station` object. If number the
            positions and offsets are retrieved from the public database.
            Otherwise the offsets will be determined with the available data.
        :param overwrite: if True overwrite existing reconstruction table.
        :param progress: if True show a progressbar while reconstructing.
        :param verbose: if True be verbose about station metadata usage.
        :param destination: alternative name for reconstruction table.

        """
        super().__init__(
            source_data, source_group, station, overwrite, progress, verbose,
            destination, force_fresh, force_stale)
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
        try:
            self.reconstructions._v_attrs.station = self.station
        except tables.HDF5ExtError:
            warnings.warn('Unable to store station object, to large for HDF.')


class ReconstructSimulatedEvents(ReconstructESDEvents):

    """Reconstruct simulated events from single stations

    Simulated events use simulated meta-data (e.g. timing offsets)
    which are stored as a :class:`~sapphire.clusters.BaseCluster` object
    The object is stored as an node attribute of '/coincidences' in the
    HDF5 file. This class will try to read that object and use it's meta-data
    in reconstructions.

    The station number must match the station number in the stored object.

    Example usage::

        >>> import tables
        >>> from sapphire import ReconstructESDEvents

        >>> data = tables.open_file('simulation.h5', 'a')
        >>> station_path = '/cluster_simulations/station_506'
        >>> rec = ReconstructESDEvents(data, station_path, 506, overwrite=True)
        >>> rec.reconstruct_and_store()
    """

    def _get_or_create_station_object(self, station):
        """
        Returns a :class:`~sapphire.clusters.Station` object that
        contains the meta data of the station. For simulated data this class
        contains the geometry of the station and detectors as well as the
        generated (random) timing offsets. This object is stored in the
        HDF5 file generated by the simulation.

        If :attr:`station` is a station number, this will try to extract
        that station from the cluster object stored in the HDF5 file. If
        :attr:`station` is a station object it will return that object.

        :param station: The station number used in the simulation or
            an :class:`~sapphire.clusters.Station` object.

        """
        if isinstance(station, Station):
            self.station = station
            self.station_number = None
            if self.verbose:
                print('Using object %s for metadata.' % self.station)
        else:
            self.station_number = station
            try:
                cluster = self.data.get_node_attr('/coincidences', 'cluster')
                self.station = cluster.get_station(station)
                if self.station is None:
                    raise RuntimeError('Station %d not found in cluster'
                                       ' object.' % self.station_number)
                if self.verbose:
                    print('Read object %s from datafile.' % self.station)
            except (tables.NoSuchNodeError, AttributeError):
                raise RuntimeError('Unable to read cluster object from HDF')


class ReconstructESDCoincidences:

    """Reconstruct coincidences, e.g. event between multiple stations

    Example usage::

        >>> import tables
        >>> from sapphire import ReconstructESDCoincidences

        >>> data = tables.open_file('2014_1_1.h5', 'a')
        >>> rec = ReconstructESDCoincidences(data, overwrite=True)
        >>> rec.reconstruct_and_store()

    """

    def __init__(self, data, coincidences_group='/coincidences',
                 overwrite=False, progress=True, verbose=False,
                 destination='reconstructions', cluster=None,
                 force_fresh=False, force_stale=False):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param coincidences_group: the destination group.
        :param overwrite: if True overwrite existing reconstruction table.
        :param progress: if True show a progressbar while reconstructing.
        :param verbose: if True be verbose about station metadata usage.
        :param destination: alternative name for reconstruction table.
        :param cluster: a Cluster object to use for the reconstructions.

        """
        self.data = data
        self.coincidences_group = data.get_node(coincidences_group)
        self.coincidences = self.coincidences_group.coincidences
        self.overwrite = overwrite
        self.progress = progress
        self.verbose = verbose
        self.destination = destination
        self.force_fresh = force_fresh
        self.force_stale = force_stale
        self.offsets = {}

        self.cq = CoincidenceQuery(data, self.coincidences_group)
        self.cluster = self._get_or_create_cluster_object(cluster)

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

        :param station_numbers: list of stations to use for reconstructions.

        """
        if len(self.core_x) and len(self.core_y):
            initials = ({'core_x': x, 'core_y': y}
                        for x, y in zip(self.core_x, self.core_y))
        else:
            initials = []
        coincidences = pbar(self.cq.all_coincidences(iterator=True),
                            length=self.coincidences.nrows, show=self.progress)
        angles = self.direction.reconstruct_coincidences(
            self.cq.all_events(coincidences, n=0), station_numbers,
            self.offsets, progress=False, initials=initials)
        self.theta, self.phi, self.station_numbers = angles

    def reconstruct_cores(self, station_numbers=None):
        """Reconstruct core for all events

        :param station_numbers: list of stations to use for reconstructions.

        """
        if len(self.theta) and len(self.phi):
            initials = ({'theta': theta, 'phi': phi}
                        for theta, phi in zip(self.theta, self.phi))
        else:
            initials = []
        coincidences = pbar(self.cq.all_coincidences(iterator=True),
                            length=self.coincidences.nrows, show=self.progress)
        cores = self.core.reconstruct_coincidences(
            self.cq.all_events(coincidences, n=0), station_numbers,
            progress=False, initials=initials)
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
        try:
            self.reconstructions._v_attrs.cluster = self.cluster
        except tables.HDF5ExtError:
            warnings.warn('Unable to store cluster object, to large for HDF.')

    def get_station_timing_offsets(self):
        """Construct a dict of :class:`~sapphire.api.Station` objects

        Try to extract offsets from provided cluster objects into a dictionary,
        to be used by the reconstructions.
        If the cluster is not available create a :class:`~sapphire.api.Station`
        object for each station in the cluster.

        """
        try:
            self.offsets = {station.number: [station.gps_offset + d.offset
                                             for d in station.detectors]
                            for station in self.cluster.stations}
            if self.verbose:
                print('Using timing offsets from cluster object.')
        except AttributeError:
            self.offsets = {station.number:
                            api.Station(station.number,
                                        force_fresh=self.force_fresh,
                                        force_stale=self.force_stale)
                            for station in self.cluster.stations}
            if self.verbose:
                print('Using timing offsets from public database.')

    def store_reconstructions(self):
        """Loop over list of reconstructed data and store results

        Unsuccessful reconstructions are also stored but with the NumPy
        NaN as reconstructed value.

        """
        for coincidence, x, y, theta, phi, station_numbers in zip_longest(
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

    def _get_or_create_cluster_object(self, cluster):
        """Create cluster object from public database"""

        if cluster is None:
            s_active = self._get_active_stations()
            cluster = HiSPARCStations(s_active,
                                      force_fresh=self.force_fresh,
                                      force_stale=self.force_stale)
            if self.verbose:
                print('Constructed cluster %s from public database.'
                      % self.cluster)
        else:
            # TODO: check cluster object isinstance
            if self.verbose:
                print('Using cluster %s for metadata.' % self.cluster)
        return cluster

    def _get_active_stations(self):
        """Return station numbers with non-empty event table in datafile"""

        active_stations = []

        for s_path in self.coincidences_group.s_index:
            try:
                station_event_table = self.data.get_node(s_path.decode() +
                                                         '/events')
            except tables.NoSuchNodeError:
                continue
            if not station_event_table.nrows:
                continue
            active_stations.append(int(s_path.split(b'station_')[-1]))

        return active_stations


class ReconstructESDCoincidencesFromSource(ReconstructESDCoincidences):

    def __init__(self, source_data, dest_data, source_group, dest_group,
                 overwrite=False, progress=True, verbose=False,
                 destination='reconstructions', cluster=None,
                 force_fresh=False, force_stale=False):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param station_group: the group containing the event table,
            the results will also be stored in this group.
        :param station: either a station number or
            :class:`sapphire.clusters.Station` object. If number the
            positions and offsets are retrieved from the public database.
            Otherwise the offsets will be determined with the available data.
        :param overwrite: if True overwrite existing reconstruction table.
        :param progress: if True show a progressbar while reconstructing.
        :param verbose: if True be verbose about station metadata usage.
        :param destination: alternative name for reconstruction table.

        """
        super().__init__(
            source_data, source_group, overwrite, progress, verbose,
            destination, cluster, force_fresh, force_stale)
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
        try:
            self.reconstructions._v_attrs.cluster = self.cluster
        except tables.HDF5ExtError:
            warnings.warn('Unable to store cluster object, to large for HDF.')


class ReconstructSimulatedCoincidences(ReconstructESDCoincidences):

    """Reconstruct simulated coincidences.

    Simulated coincidences use simulated meta-data (e.g. timing offsets)
    which are stored as a :class:`~sapphire.clusters.BaseCluster` object
    The object is stored as an node attribute of '/coincidences' in the
    HDF5 file. This class will try to read that object and use it's meta-data
    in reconstructions.

    Example usage::

        >>> import tables
        >>> from sapphire import ReconstructSimulatedCoincidences

        >>> data = tables.open_file('simulated.h5', 'a')
        >>> rec = ReconstructSimulatedCoincidences(data, overwrite=True)
        >>> rec.reconstruct_and_store()

    """

    def _get_or_create_cluster_object(self, cluster):
        """
        Returns a :class:`~sapphire.clusters.BaseCluster` object that
        contains the meta data of the cluster. For simulated data this class
        contains the geometry of the stations and detectors as well as the
        generated (random) timing offsets. This object stored in the HDF5 file
        generated by the simulation.

        If :attr:`cluster` is `None`, this will try to extract the cluster
        object from the HDF5 file. If :attr:`cluster` is a cluster
        object it will return that object.

        :param cluster: a :class:`~sapphire.clusters.BaseCluster` object or
            `None`.

        """
        if cluster is None:
            try:
                cluster = self.data.get_node_attr(self.coincidences_group,
                                                  'cluster')
                if self.verbose:
                    print('Read cluster %s from datafile.' % self.cluster)
            except (tables.NoSuchNodeError, AttributeError):
                raise RuntimeError('Unable to read cluster object from HDF')
        else:
            # TODO: check cluster object
            if self.verbose:
                print('Using cluster %s for metadata.' % self.cluster)
        return cluster
