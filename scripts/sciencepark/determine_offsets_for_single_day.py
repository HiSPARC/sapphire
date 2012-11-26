from math import sqrt
import datetime
import operator
import os.path
import itertools

import tables
import numpy as np
from numpy import arctan2, cos, sin, arcsin, isnan, pi
import progressbar as pb

from hisparc.publicdb import download_data
from sapphire.analysis.process_events import ProcessEventsWithLINT
from sapphire import storage, clusters


class Master:
    stations = range(501, 503)
    datetimerange = (datetime.datetime(2012, 3, 2),
                     datetime.datetime(2012, 3, 3))

    def __init__(self, data_path):
        self.data = tables.openFile(data_path, 'a')

        self.station_groups = ['/s%d' % u for u in self.stations]
        self.cluster = clusters.ScienceParkCluster(self.stations)

        self.trig_threshold = .5

    def main(self):
        self.download_data()
        self.clean_data()
        self.process_events()

    def download_data(self):
        start, end = self.datetimerange

        for station, group_path in zip(self.stations, self.station_groups):
            if not group_path in self.data:
                print "Downloading data for station", station
                download_data(self.data, group_path, station,
                              start, end, get_blobs=True)

    def clean_data(self):
        for group in self.station_groups:
            group = self.data.getNode(group)
            attrs = group._v_attrs
            if not 'is_clean' in attrs or not attrs.is_clean:
                self.clean_events_in_group(group)
                attrs.is_clean = True

    def clean_events_in_group(self, group):
        events = group.events

        timestamps = [u for u in enumerate(events.col('ext_timestamp'))]
        timestamps.sort(key=operator.itemgetter(1))

        prev = 0
        unique_ids = []
        for row_id, timestamp in timestamps:
            if timestamp != prev:
                unique_ids.append(row_id)
            prev = timestamp

        tmptable = self.data.createTable(group, 't__events',
                                         description=events.description)
        rows = events.readCoordinates(unique_ids)
        tmptable.append(rows)
        tmptable.flush()
        self.data.renameNode(tmptable, 'events', overwrite=True)

    def process_events(self):
        attrs = self.data.root._v_attrs
        for station_id, station_group in enumerate(self.station_groups):
            process = ProcessEventsWithLINT(self.data, station_group)
            if 'is_processed' not in attrs or not attrs.is_processed:
                process.process_and_store_results()
            offsets = process.determine_detector_timing_offsets()
            print "Offsets for station %d: %s" % (station_id, offsets)

        attrs.is_processed = True


if __name__ == '__main__':
#    np.seterr(divide='ignore', invalid='ignore')

    master = Master('master-singleday.h5')
    master.main()
