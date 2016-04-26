""" Determine time differences between coincident events

Determine time delta between coincidence events from station pairs.

"""
import tables
from numpy import isnan

from .coincidence_queries import CoincidenceQuery
from .event_utils import station_arrival_time
from ..storage import TimeDelta


class ProcessTimeDeltas(object):

    """Process HiSPARC event coincidences to obtain time deltas.

    Use this to determine arrival time differences between station pairs which
    have coincident events.

    """

    def __init__(self, data, coincidence_group='/coincidences', progress=True):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param coincidence_group: path to the coincidences group.
        :param progress: show progressbar.

        """
        self.data = data
        self.cq = CoincidenceQuery(self.data, coincidence_group)
        self.progress = progress

    def determine_and_store_time_deltas(self):
        self.determine_time_differences()
        self.store_time_deltas()

    def determine_time_deltas(self, coin_events, ref_station, station,
                                   ref_detector_offsets=None,
                                   detector_offsets=None):
        """Determine the arrival time differences between two stations.

        :param coin_events: coincidence events from a CoincidenceQuery.
        :param ref_station,station: station numbers.
        :param ref_detector_offsets,detector_offsets: detector timing offset
            list. If None the station numbers are used to get the offsets from
            the API.
        :return: extended timestamp of the first event and time difference,
                 t - t_ref. Not corrected for altitude differences.

        """
        dt = []
        ets = []
        previous_ets = 0

        if ref_detector_offsets is None:
            Station(ref_station).detector_timing_offsets
        if detector_offsets is None:
            Station(station).detector_timing_offsets

        for events in coin_events:
            ref_ets = events[0][1]['ext_timestamp']
            # Filter coincidence which is subset of previous coincidence
            if previous_ets == ref_ets:
                continue
            else:
                previous_ets = ref_ets
            # Filter for possibility of same station twice in coincidence
            if len(events) is not 2:
                continue
            if events[0][0] == ref_station:
                ref_id = 0
                id = 1
            else:
                ref_id = 1
                id = 0
            ref_t = station_arrival_time(events[ref_id][1], ref_ets, [0, 1, 2, 3],
                                         ref_detector_offsets)
            t = station_arrival_time(events[id][1], ref_ets, [0, 1, 2, 3],
                                     detector_offsets)
            if isnan(t) or isnan(ref_t):
                continue
            dt.append(t - ref_t)
            ets.append(ref_ets)
        return ets, dt

    def store_time_deltas(self, data, ref_station, station, ext_timestamps,
                          time_deltas):
        """Store determined dt values"""

        table_path = '/time_deltas/station_%d/station_%d' % (ref_station, station)
        try:
            dt_table = data.get_node(table_path, 'time_deltas')
            dt_table.remove()
        except tables.NoSuchNodeError:
            pass
        delta_data = [(ets, int(ets) / int(1e9), int(ets) % int(1e9), time_delta)
                      for ets, time_delta in zip(ext_timestamps, time_deltas)]
        table = data.create_table(table_path, 'time_deltas', TimeDelta,
                                  createparents=True, expectedrows=len(delta_data))
        table.append(delta_data)
        table.flush()
