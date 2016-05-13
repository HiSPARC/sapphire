""" Determine time differences between coincident events

    Determine time delta between coincidence events from station pairs.

    Example usage::

        import datetime

        import tables

        from sapphire import download_coincidences
        from sapphire import ProcessTimeDeltas

        START = datetime.datetime(2015, 2, 1)
        END = datetime.datetime(2015, 2, 5)

        if __name__ == '__main__':
            with tables.open_file('data.h5', 'w') as data:
                download_coincidences(data, start=START, end=END)
                td = ProcessTimeDeltas(data)
                td.determine_and_store_time_deltas()

"""
import re
from itertools import combinations

import tables
from numpy import isnan

from ..utils import pbar
from ..api import Station
from ..storage import TimeDelta
from .coincidence_queries import CoincidenceQuery
from .event_utils import station_arrival_time


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
        """Find station pairs, determine time deltas, and store the results."""

        self.find_station_pairs()
        self.get_detector_offsets()
        self.determine_and_store_time_deltas_for_pairs()

    def determine_and_store_time_deltas_for_pairs(self):
        """Determine time deltas for all pairs and store the results."""

        for pair in pbar(self.pairs, show=self.progress):
            ets, dt = self.determine_time_deltas_for_pair(*pair)
            if len(ets):
                self.store_time_deltas(ets, dt, pair)

    def find_station_pairs(self):
        """Find all unique station pairs which are in a coincidence together

        Assumes the stations in the s_index are sorted by station number.

        """
        s_index = self.cq.s_index
        re_number = re.compile('[0-9]+$')
        s_numbers = [int(re_number.search(s_path).group())
                     for s_path in s_index]

        c_index = self.cq.c_index
        self.pairs = {(s_numbers[s1], s_numbers[s2])
                      for c_idx in c_index
                      for s1, s2 in combinations(sorted(c_idx[:, 0]), 2)}

    def get_detector_offsets(self):
        """Retrieve the API detector_timing_offset method for all pairs

        The detector_timing_offset methods accept a single timestamp as
        argument, and return the detector offsets for this timestamp.

        """
        station_numbers = {station for pair in self.pairs for station in pair}
        self.detector_timing_offsets = {sn: Station(sn).detector_timing_offset
                                        for sn in station_numbers}

    def determine_time_deltas_for_pair(self, ref_station, station):
        """Determine the arrival time differences between two stations.

        :param ref_station,station: station numbers.
        :return: extended timestamp of the first event and time difference,
                 t - t_ref. Not corrected for altitude differences.

        """
        dt = []
        ets = []
        previous_ets = 0

        coincidences = self.cq.all([ref_station, station], iterator=True)
        coin_events = self.cq.events_from_stations(coincidences,
                                                   [ref_station, station])

        ref_offsets = self.detector_timing_offsets[ref_station]
        offsets = self.detector_timing_offsets[station]

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

            ref_event = events[ref_id][1]
            ref_detector_offsets = ref_offsets(ref_event['timestamp'])
            event = events[id][1]
            detector_offsets = offsets(event['timestamp'])

            ref_t = station_arrival_time(ref_event, ref_ets, [0, 1, 2, 3],
                                         ref_detector_offsets)
            t = station_arrival_time(event, ref_ets, [0, 1, 2, 3],
                                     detector_offsets)
            if isnan(t) or isnan(ref_t):
                continue
            dt.append(t - ref_t)
            ets.append(ref_ets)
        return ets, dt

    def store_time_deltas(self, ext_timestamps, time_deltas, pair):
        """Store determined dt values"""

        table_path = '/time_deltas/station_%d/station_%d' % pair
        try:
            dt_table = self.data.get_node(table_path, 'time_deltas')
            dt_table.remove()
        except tables.NoSuchNodeError:
            pass
        delta_data = [(ets, int(ets) / int(1e9), int(ets) % int(1e9),
                       time_delta)
                      for ets, time_delta in zip(ext_timestamps, time_deltas)]
        table = self.data.create_table(table_path, 'time_deltas', TimeDelta,
                                       createparents=True,
                                       expectedrows=len(delta_data))
        table.append(delta_data)
        table.flush()
