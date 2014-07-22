import itertools

import progressbar
from numpy import nan, isnan, arange, histogram
from scipy.optimize import curve_fit
from scipy.stats import norm
import tables

from ..storage import ReconstructedEvent, ReconstructedCoincidence
from ..clusters import HiSPARCStations, ScienceParkCluster
from .direction_reconstruction import (DirectEventReconstruction,
                                       FitEventReconstruction,
                                       DirectCoincidenceReconstruction,
                                       FitCoincidenceReconstruction)
from .coincidence_queries import CoincidenceQuery


class ReconstructESDEvents(object):

    """Reconstruct events from single stations

    Example usage::

        import tables
        from sapphire.analysis.reconstructions import ReconstructESDEvents

        data = tables.open_file('2014_1_1.h5', 'a')
        station_path = '/hisparc/cluster_amsterdam/station_506'
        dirrec = ReconstructESDEvents(data, station_path, 506, overwrite=True)
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
        self.offsets = [0., 0., 0., 0.]

        self.station = (HiSPARCStations([station_number])
                        .get_station(station_number))

        self.direct = DirectEventReconstruction(self.station)
        self.fit = FitEventReconstruction(self.station)

    def reconstruct_and_store(self):
        """Shorthand function to reconstruct event and store the results"""

        self.prepare_output()
        self.determine_detector_timing_offsets()
        self.store_offsets()
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
            theta, phi = self.direct.reconstruct_event(event, detector_ids,
                                                       self.offsets)
        elif len(detector_ids) == 4:
            theta, phi = self.fit.reconstruct_event(event, self.offsets)
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

    def determine_detector_timing_offsets(self):
        """Determine the offsets between the station detectors.

        ADL: Currently assumes detector 1 is a good reference.
        But this is not always the best choice. Perhaps it should be
        determined using more data (more than one day) to be more
        accurate. Also assumes the detectors are at the same altitude.

        """
        gauss = lambda x, N, m, s: N * norm.pdf(x, m, s)
        bins = arange(-100 + 1.25, 100, 2.5)

        t2 = self.events.col('t2')

        offsets = []
        for timings in 't1', 't3', 't4':
            timings = self.events.col(timings)
            dt = (timings - t2).compress((t2 >= 0) & (timings >= 0))
            y, bins = histogram(dt, bins=bins)
            x = (bins[:-1] + bins[1:]) / 2
            try:
                popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0., 10.))
                offsets.append(popt[1])
            except RuntimeError:
                offsets.append(0.)

        self.offsets = offsets[0:1] + [0.] + offsets[1:]

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
            row['d%d' % (id + 1)] = True
        row.append()


class ReconstructESDCoincidences(object):

    """Reconstruct coincidences, e.g. event between multiple stations

    Example usage::

        import tables
        from sapphire.analysis.reconstructions import ReconstructESDCoincidences

        data = tables.open_file('2014_1_1.h5', 'a')
        dirrec = ReconstructESDCoincidences(data, overwrite=True)
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
        self.offsets = {}

        self.cq = CoincidenceQuery(data, self.coincidences_group)
        # Get latest position data
        s_numbers = [station.number for station in
                     self.coincidences_group._f_getattr('cluster').stations]
        self.cluster = HiSPARCStations(s_numbers)

        self.direct = DirectCoincidenceReconstruction(self.cluster)
        self.fit = FitCoincidenceReconstruction(self.cluster)

    def reconstruct_and_store(self):
        """Shorthand function to reconstruct coincidences and store results"""

        self.prepare_output()
        self.determine_station_timing_offsets()
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
            theta, phi = self.direct.reconstruct_coincidence(coincidence,
                                                             self.offsets)
        elif len(coincidence) >= 4:
            theta, phi = self.fit.reconstruct_coincidence(coincidence,
                                                          self.offsets)
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

    def determine_station_timing_offsets(self):
        """Determine the offsets between the stations."""

        for s_path in self.coincidences_group.s_index:
            station_number = int(s_path.split('station_')[-1])
            station_group = self.data.get_node(s_path)
            offsets = self.determine_detector_timing_offsets(station_group)
            self.offsets[station_number] = offsets

    def determine_detector_timing_offsets(self, station_group):
        """Determine the offsets between the station detectors.

        ADL: Currently assumes detector 1 is a good reference.
        But this is not always the best choice. Perhaps it should be
        determined using more data (more than one day) to be more
        accurate. Also assumes the detectors are at the same altitude.

        """
        gauss = lambda x, N, m, s: N * norm.pdf(x, m, s)
        bins = arange(-100 + 1.25, 100, 2.5)

        t2 = station_group.events.col('t2')

        offsets = []
        for timings in 't1', 't3', 't4':
            timings = station_group.events.col(timings)
            dt = (timings - t2).compress((t2 >= 0) & (timings >= 0))
            y, bins = histogram(dt, bins=bins)
            x = (bins[:-1] + bins[1:]) / 2
            try:
                popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0., 10.))
                offsets.append(popt[1])
            except RuntimeError:
                offsets.append(0.)
        return offsets[0:1] + [0.] + offsets[1:]

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

        row['reference_x'] = coincidence['x']
        row['reference_y'] = coincidence['y']
        row['reference_zenith'] = coincidence['zenith']
        row['reference_azimuth'] = coincidence['azimuth']
        row['reference_size'] = coincidence['size']
        row['reference_energy'] = coincidence['energy']

        for number in station_numbers:
            row['s%d' % number] = True

        row.append()


def pbar(iterator):
    """Get a new progressbar with our default widgets"""

    pb = progressbar.ProgressBar(widgets=[progressbar.ETA(), progressbar.Bar(),
                                          progressbar.Percentage()])
    return pb(iterator)
