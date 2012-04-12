import datetime
import operator

import tables
import numpy as np
import progressbar as pb

from hisparc.publicdb import download_data
from hisparc.analysis import coincidences
from sapphire.analysis.process_events import ProcessIndexedEvents
from sapphire import storage, clusters
import transformations


class ScienceParkCluster(clusters.BaseCluster):
    # website 1 maart 2012
#    gps_coordinates = {501: (52.3559126101, 4.95074843056,
#                               58.6631137794),
#                       503: (52.3562664259, 4.95294989286,
#                               48.7995043574),
#                       506: (52.3571869092, 4.95198066585,
#                               45.4739832897),
#                      }

    # Niels mailtje 6 nov 2011
#    gps_coordinates = {501: (52.3559179545,4.95114534876,
#                               58.6631137794),
#                       503: (52.3562664259,4.95294989286,
#                               48.7995043574),
#                       506: (52.3571787512,4.95198605591,
#                               45.4739832897),
#                      }

    # 1 day self-survey (8 april 2011) + 506 (Niels, pos from site on
    # 2 dec, 2011)
    gps_coordinates = {501: (52.355924173294305, 4.951144021644267,
                             56.102345941588283),
                       502: (52.355293344895919, 4.9501047083812697,
                             55.954367009922862),
                       503: (52.356254735127557, 4.9529437445598328,
                             51.582641703076661),
                       504: (52.357178777910278, 4.9543838852175561,
                             54.622688433155417),
                       505: (52.357251580629246, 4.9484007564706891,
                             47.730995402671397),
                       506: (52.3571787512,4.95198605591,
                             43.8700314863),
                      }

    station_rotations = {501: 135, 503: 45, 506: 267}


    def __init__(self, stations):
        super(ScienceParkCluster, self).__init__()

        reference = self.gps_coordinates[501]
        transformation = \
            transformations.FromWGS84ToENUTransformation(reference)

        for station in stations:
            easting, northing, up = \
                transformation.transform(self.gps_coordinates[station])
            alpha = self.station_rotations[station] / 180. * np.pi
            self._add_station((easting, northing), alpha)


class Master:
    stations = [501, 503, 506]
    datetimerange = (datetime.datetime(2012, 3, 1),
                     datetime.datetime(2012, 3, 2))

    def __init__(self, data_path):
        self.data = tables.openFile('master.h5', 'a')

        self.station_groups = ['/s%d' % u for u in self.stations]
        self.cluster = ScienceParkCluster(self.stations)

        self.trig_threshold = .5

    def main(self):
        self.download_data()
        self.clean_data()

        self.search_coincidences()
        self.process_events_from_c_index()
        #self.data.removeNode('/coincidences', recursive=True)
        self.store_coincidences()

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

    def search_coincidences(self):
        if '/c_index' not in self.data and '/timestamps' not in self.data:
            c_index, timestamps = \
                coincidences.search_coincidences(self.data,
                                                 self.station_groups)
            timestamps = np.array(timestamps, dtype=np.uint64)
            self.data.createArray('/', 'timestamps', timestamps)
            self.data.createVLArray('/', 'c_index', tables.UInt32Atom())
            for coincidence in c_index:
                self.data.root.c_index.append(coincidence)

    def process_events_from_c_index(self):
        attrs = self.data.root._v_attrs
        if 'is_processed' not in attrs or not attrs.is_processed:
            c_index = self.data.root.c_index.read()
            timestamps = self.data.root.timestamps.read()

            selected_timestamps = []
            for coincidence in c_index:
                for event in coincidence:
                    selected_timestamps.append(timestamps[event])
            full_index = np.array(selected_timestamps)

            for station_id, station_group in enumerate(self.station_groups):
                selected = full_index.compress(full_index[:, 1] == station_id,
                                               axis=0)
                index = selected[:, 2]

                process = ProcessIndexedEvents(self.data, station_group,
                                               index)
                process.process_and_store_results()

            attrs.is_processed = True

    def store_coincidences(self):
        if '/coincidences' not in self.data:
            group = self.data.createGroup('/', 'coincidences')
            self.c_index = []
            self.coincidences = self.data.createTable(group,
                                                      'coincidences',
                                                      storage.Coincidence)
            self.observables = self.data.createTable(group, 'observables',
                                            storage.EventObservables)

            progress = pb.ProgressBar(widgets=[pb.ETA(), pb.Bar(),
                                               pb.Percentage()])
            for coincidence in progress(self.data.root.c_index):
                self.store_coincidence(coincidence)

            c_index = self.data.createVLArray(group, 'c_index',
                                              tables.UInt32Col())
            for coincidence in self.c_index:
                c_index.append(coincidence)
            c_index.flush()
            self.c_index = c_index

    def store_coincidence(self, coincidence):
        row = self.coincidences.row
        coincidence_id = len(self.coincidences)
        row['id'] = coincidence_id
        row['N'] = len(coincidence)

        observables_idx = []
        timestamps = []
        for index in coincidence:
            event_desc = self.data.root.timestamps[index]
            station_id = event_desc[1]
            event_index = event_desc[2]

            group = self.data.getNode(self.station_groups[station_id])
            event = group.events[event_index]
            idx = self.store_event_in_observables(event, coincidence_id,
                                                  station_id)
            observables_idx.append(idx)
            timestamps.append((event['ext_timestamp'], event['timestamp'],
                               event['nanoseconds']))

        first_timestamp = sorted(timestamps)[0]
        row['ext_timestamp'], row['timestamp'], row['nanoseconds'] = \
            first_timestamp
        row.append()
        self.c_index.append(observables_idx)
        self.coincidences.flush()

    def store_event_in_observables(self, event, coincidence_id, station_id):
        row = self.observables.row
        event_id = len(self.observables)
        row['id'] = event_id

        row['station_id'] = station_id
        for key in ('timestamp', 'nanoseconds', 'ext_timestamp',
                    'n1', 'n2', 'n3', 'n4', 't1', 't2', 't3', 't4'):
            row[key] = event[key]

        signals = [event[key] for key in 'n1', 'n2', 'n3', 'n4']
        N = sum([1 if u > self.trig_threshold else 0 for u in signals])
        row['N'] = N

        row.append()
        self.observables.flush()
        return event_id


if __name__ == '__main__':
    master = Master('data.h5')
    master.main()
