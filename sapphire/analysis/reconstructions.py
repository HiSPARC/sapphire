import itertools

import progressbar
from numpy import nan, isnan

from ..storage import ReconstructedEvent, ReconstructedCoincidence
from ..clusters import HiSPARCStations, ScienceParkCluster
from .direction_reconstruction import (DirectEventReconstruction,
                                       FitEventReconstruction)


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
        self.station_group = data.getNode(station_group)
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
        row['zenith'] = theta
        row['azimuth'] = phi
        for id in detector_ids:
            row['d%d' % id] = True
        row.append()


def pbar(iterator):
    """Get a new progressbar with our default widgets"""

    pb = progressbar.ProgressBar(widgets=[progressbar.ETA(), progressbar.Bar(),
                                          progressbar.Percentage()])
    return pb(iterator)
