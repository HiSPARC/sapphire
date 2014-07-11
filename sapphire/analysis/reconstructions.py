import itertools

import progressbar
from numpy import nan, isnan
import tables

from ..storage import ReconstructedEvent, ReconstructedCoincidence
from ..clusters import HiSPARCStations, ScienceParkCluster
from .direction_reconstruction import (DirectEventReconstruction,
                                       FitEventReconstruction,
                                       DirectClusterReconstruction,
                                       FitClusterReconstruction)
from .coincidence_queries import CoincidenceQuery


class ReconstructESDEvents(object):

    """Reconstruct events from single stations

    Example usage::

        import tables
        from sapphire.analysis.reconstructions import ReconstructESDEvents

        data = tables.open_file('2014/1/2014_1_2.h5', 'a')
        station_path = '/hisparc/cluster_amsterdam/station_506'
        dirrec = ReconstructESDEvents(data, station_path, 506)
        dirrec.reconstruct_and_store()

    """

    def __init__(self, data, station_group, station_number,
                 overwrite=False, progress=True):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param station_group: the destination group.
        :param station_number: station identifier.
        :param overwrite: if True, overwrite existing reconstruction table.
        :param progress: if True, show a progressbar while reconstructing.

        """
        self.data = data
        self.station_group = data.get_node(station_group)
        self.events = self.station_group.events
        self.overwrite = overwrite
        self.progress = progress

        self.station = (HiSPARCStations([station_number])
                        .get_station(station_number))

        self.direct = DirectEventReconstruction(self.station)
        self.fit = FitEventReconstruction(self.station)

    def reconstruct_and_store(self):
        """Shorthand function to reconstruct event and store the results"""

        self.prepare_output()
        self.reconstruct_directions()
        self.store_reconstructions()

    def reconstruct_directions(self):
        """Reconstruct all events

        Reconstruct each event in the events tables.

        """
        events = pbar(self.events) if self.progress else self.events
        angles = [self._reconstruct_direction(e) for e in events]
        self.theta, self.phi, self.detector_ids = zip(*angles)

    def _reconstruct_direction(self, event):
        """Reconstruct an event

        Use direct algorithm if three detectors have an arrival time,
        use fit algorithm in case of four and return (nan, nan) otherwise.

        """
        detector_ids = [id for id in range(4)
                        if event['t%d' % (id + 1)] not in [-1, -999]]
        if len(detector_ids) == 3:
            theta, phi = self.direct.reconstruct_event(event, detector_ids)
        elif len(detector_ids) == 4:
            theta, phi = self.fit.reconstruct_event(event)
        else:
            theta, phi = (nan, nan)
        return theta, phi, detector_ids

    def prepare_output(self):
        """Prepare output table"""

        if 'reconstructions' in self.station_group:
            if self.overwrite:
                self.data.remove_node(self.station_group.reconstructions,
                                      recursive=True)
            else:
                raise RuntimeError("Reconstructions table already exists for "
                                   "%s, and overwrite is False" %
                                   self.station_group)
        self.reconstructions = self.data.create_table(
            self.station_group, 'reconstructions', ReconstructedEvent)
        self.reconstructions._v_attrs.station = self.station

    def store_reconstructions(self):
        """Loop over list of reconstructed data and store results

        Only writes rows if reconstruction was possible and successful.

        ADL: Perhaps we should always store reconstructions, and use
        error values in case it failed. However, the usual -999 might be
        a real value (though unlikely to be exactly -999) in case of
        core position reconstruction.

        """
        for event, theta, phi, detector_ids in itertools.izip(
                self.events, self.theta, self.phi, self.detector_ids):
            if not isnan(theta) and not isnan(phi):
                self._store_reconstruction(event, theta, phi, detector_ids)
        self.reconstructions.flush()

    def _store_reconstruction(self, event, theta, phi, detector_ids):
        """Store single reconstruction"""

        row = self.reconstructions.row
        row['id'] = event['event_id']
        row['ext_timestamp'] = event['ext_timestamp']
        row['min_n'] = min([event['n%d' % (id + 1)] for id in detector_ids])
        row['zenith'] = theta
        row['azimuth'] = phi
        for id in detector_ids:
            row['d%d' % id] = True
        row.append()


class ReconstructESDCoincidences(object):

    """Reconstruct coincidences, e.g. event between multiple stations

    Example usage::

        import tables
        from sapphire.analysis.reconstructions import ReconstructESDCoincidences

        data = tables.open_file('2014/1/2014_1_2.h5', 'a')
        dirrec = ReconstructESDCoincidences(data)
        dirrec.reconstruct_and_store()

    """

    def __init__(self, data, coincidences_group='/coincidences',
                 overwrite=False, progress=True):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param coincidences_group: the destination group.
        :param overwrite: if True, overwrite existing reconstruction table.
        :param progress: if True, show a progressbar while reconstructing.

        """
        self.data = data
        self.coincidences_group = data.get_node(coincidences_group)
        self.coincidences = self.coincidences_group.coincidences
        self.overwrite = overwrite
        self.progress = progress
        self.cq = CoincidenceQuery(data, self.coincidences_group)
        # Get latest position data
        s_numbers = [station.number for station in
                     self.coincidences_group._f_getattr('cluster').stations]
        self.cluster = HiSPARCStations(s_numbers)

        self.direct = DirectClusterReconstruction(self.cluster)
        self.fit = FitClusterReconstruction(self.cluster)

    def reconstruct_and_store(self):
        """Shorthand function to reconstruct coincidences and store results"""

        self.prepare_output()
        self.reconstruct_directions()
        self.store_reconstructions()

    def reconstruct_directions(self):
        """Reconstruct all coincidences

        Reconstruct each coincidence in the coincidences tables.

        """
        coincidences = self.cq.all_coincidences()
        coincidence_events = self.cq.all_events(coincidences)
        if self.progress:
            coincidence_events = pbar(coincidence_events)
        angles = [self._reconstruct_direction(c) for c in coincidence_events]
        self.theta, self.phi, self.station_numbers = zip(*angles)

    def _reconstruct_direction(self, coincidence):
        """Reconstruct a coincidence

        Use direct algorithm if three stations are in coincidence,
        use fit algorithm in case of four or more,
        return (nan, nan) otherwise.

        """
        station_numbers = [c[0] for c in coincidence]
        if len(coincidence) == 3:
            theta, phi = self.direct.reconstruct_coincidence(coincidence)
        elif len(coincidence) >= 4:
            theta, phi = self.fit.reconstruct_coincidence(coincidence)
        else:
            theta, phi = (nan, nan)
        return theta, phi, station_numbers

    def prepare_output(self):
        """Prepare output table"""

        if 'reconstructions' in self.coincidences_group:
            if self.overwrite:
                self.data.remove_node(self.coincidences_group.reconstructions,
                                      recursive=True)
            else:
                raise RuntimeError("Reconstructions table already exists for "
                                   "%s, and overwrite is False" %
                                   self.coincidences_group)

        s_columns = {'s%d' % station.number: tables.BoolCol(pos=p)
                     for p, station in enumerate(self.cluster.stations, 26)}
        description = ReconstructedCoincidence
        description.columns.update(s_columns)
        self.reconstructions = self.data.create_table(
            self.coincidences_group, 'reconstructions', description)
        self.reconstructions._v_attrs.cluster = self.cluster

    def store_reconstructions(self):
        """Loop over list of reconstructed data and store results

        Only writes rows if reconstruction was possible and successful.

        ADL: Perhaps we should always store reconstructions, and use
        error values in case it failed. However, the usual -999 might be
        a real value (though unlikely to be exactly -999) in case of
        core position reconstruction.

        """
        for coincidence, theta, phi, station_numbers in itertools.izip(
                self.coincidences, self.theta, self.phi, self.station_numbers):
            if not isnan(theta) and not isnan(phi):
                self._store_reconstruction(coincidence, theta, phi,
                                           station_numbers)
        self.reconstructions.flush()

    def _store_reconstruction(self, coincidence, theta, phi, station_numbers):
        """Store single reconstruction"""

        row = self.reconstructions.row
        row['id'] = coincidence['id']
        row['ext_timestamp'] = coincidence['ext_timestamp']
        row['zenith'] = theta
        row['azimuth'] = phi
        for number in station_numbers:
            row['s%d' % number] = True
        row.append()


def pbar(iterator):
    """Get a new progressbar with our default widgets"""

    pb = progressbar.ProgressBar(widgets=[progressbar.ETA(), progressbar.Bar(),
                                          progressbar.Percentage()])
    return pb(iterator)
