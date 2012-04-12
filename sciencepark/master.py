import datetime

import tables
import numpy as np
import progressbar as pb

from hisparc.publicdb import download_data
from hisparc.analysis import coincidences
from sapphire.analysis.process_events import ProcessIndexedEvents
from sapphire import storage


class Master:
    stations = [501, 503, 506]
    datetimerange = (datetime.datetime(2012, 3, 1),
                     datetime.datetime(2012, 3, 2))

    def __init__(self, data_path):
        self.data = tables.openFile('master.h5', 'a')

        self.station_groups = ['/s%d' % u for u in self.stations]

        self.trig_threshold = .5

    def main(self):
        self.download_data()
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
                                            storage.SimulationEventObservables)

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
        row.append()
        self.coincidences.flush()

        observables_idx = []
        for index in coincidence:
            event_desc = self.data.root.timestamps[index]
            station_id = event_desc[1]
            event_index = event_desc[2]

            group = self.data.getNode(self.station_groups[station_id])
            event = group.events[event_index]
            idx = self.store_event_in_observables(event, coincidence_id,
                                                  station_id)
            observables_idx.append(idx)
        self.c_index.append(observables_idx)

    def store_event_in_observables(self, event, coincidence_id, station_id):
        row = self.observables.row
        event_id = len(self.observables)
        row['id'] = event_id

        row['station_id'] = station_id
        for key in 'n1', 'n2', 'n3', 'n4', 't1', 't2', 't3', 't4':
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
