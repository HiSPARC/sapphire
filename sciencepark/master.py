from math import sqrt
import datetime
import operator
import os.path
import itertools

import tables
import numpy as np
from numpy import arctan2, cos, sin, arcsin, isnan, pi, linspace
import progressbar as pb
from scipy.optimize import curve_fit
from scipy.stats import norm

from hisparc.publicdb import download_data
from hisparc.analysis import coincidences
from sapphire.analysis.process_events import ProcessEvents, ProcessIndexedEventsWithLINT
from sapphire.analysis.direction_reconstruction import \
        DirectionReconstruction
from sapphire import storage, clusters


class Master:
    stations = range(501, 507)
    datetimerange = (datetime.datetime(2011, 5, 13),
                     datetime.datetime(2012, 3, 1))


    def __init__(self, data_path):
        self.data = tables.openFile(data_path, 'a')

        self.station_groups = ['/s%d' % u for u in self.stations]
        self.cluster = clusters.ScienceParkCluster(self.stations)

        self.trig_threshold = .5

        self.detector_offsets = []
        self.station_offsets = []

    def main(self):
        self.download_data()
        self.clean_data()

        self.search_coincidences()
        self.process_events_from_c_index()
        self.store_coincidences()

        self.determine_detector_offsets()
        self.determine_station_offsets()

        self.reconstruct_direction()

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

                process = ProcessIndexedEventsWithLINT(self.data, station_group,
                                                       index)
                process.process_and_store_results()

            attrs.is_processed = True

    def store_coincidences(self):
        if '/coincidences' not in self.data:
            group = self.data.createGroup('/', 'coincidences')
            group._v_attrs.cluster = self.cluster

            self.c_index = []
            self.coincidences = self.data.createTable(group,
                                                      'coincidences',
                                                      storage.Coincidence)
            self.observables = self.data.createTable(group, 'observables',
                                            storage.EventObservables)

            progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                               pb.ETA()])
            for coincidence in progress(self.data.root.c_index):
                self.store_coincidence(coincidence)

            c_index = self.data.createVLArray(group, 'c_index',
                                              tables.UInt32Col())
            for coincidence in self.c_index:
                c_index.append(coincidence)
            c_index.flush()
            self.c_index = c_index
        else:
            # Force new cluster geometry
            group = self.data.getNode('/', 'coincidences')
            group._v_attrs.cluster = self.cluster

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

    def reconstruct_direction(self):
        reconstruction = ClusterDirectionReconstruction(self.data,
                            self.stations, '/reconstructions',
                            detector_offsets=self.detector_offsets,
                            station_offsets=self.station_offsets,
                            overwrite=True)
        reconstruction.reconstruct_angles('/coincidences')

    def determine_detector_offsets(self):
        for station_id, station_group in enumerate(self.station_groups):
            process = ProcessEvents(self.data, station_group)
            offsets = process.determine_detector_timing_offsets()
            # FIXME: ugly hack for 501
            #if station_id == 0:
            #    offsets = [0, 0, 0, 0]
            print "Offsets for station %d: %s" % (station_id, offsets)
            self.detector_offsets.append(offsets)

    def determine_station_offsets(self):
        ref_group = '/s501'
        station_groups = list(self.station_groups)
        station_groups.remove(ref_group)

        gauss = lambda x, N, mu, sigma: N * norm.pdf(x, mu, sigma)
        bins = linspace(-1e3, 1e3, 101)

        for station_id, station_group in enumerate(station_groups):
            c_index, timestamps = coincidences.search_coincidences(
                                    self.data, [ref_group, station_group])

            dt = []
            c_index = [c for c in c_index if len(c) == 2]
            for i, j in c_index:
                stations = [timestamps[u][1] for u in [i, j]]
                t0, t1 = [int(timestamps[u][0]) for u in [i, j]]
                if stations[0] > stations[1]:
                    t0, t1 = t1, t0
                dt.append(t1 - t0)
            print ref_group, station_group, len(dt),
            y, bins = np.histogram(dt, bins=bins)
            x = (bins[:-1] + bins[1:]) / 2
            popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0, 100.))
            print popt
            self.station_offsets.append(popt[1])

        ref_idx = self.station_groups.index(ref_group)
        self.station_offsets.insert(ref_idx, 0.)


class ClusterDirectionReconstruction(DirectionReconstruction):
    reconstruction_description = {'coinc_id': tables.UInt32Col(),
                                  'N': tables.UInt8Col(),
                                  'reconstructed_theta': tables.Float32Col(),
                                  'reconstructed_phi': tables.Float32Col(),
                                  'min_n134': tables.Float32Col(),
                                 }
    reconstruction_coincidence_description = {'id': tables.UInt32Col(),
                                              'N': tables.UInt8Col(),
                                             }


    def __init__(self, datafile, stations, results_group=None,
                 min_n134=1., N=None, detector_offsets=None,
                 station_offsets=None, overwrite=False):
        self.data = datafile
        self.stations = stations

        if results_group:
            self.results_group = self._create_reconstruction_group_and_tables(results_group, overwrite)
        else:
            self.results_group = None

        self.min_n134 = min_n134
        self.N = N
        self.detector_offsets = detector_offsets
        self.station_offsets = station_offsets

    def _create_reconstruction_group_and_tables(self, results_group, overwrite):
        if results_group in self.data:
            if overwrite:
                self.data.removeNode(results_group, recursive=True)
            else:
                raise RuntimeError("Result group exists, but overwrite is False")

        head, tail = os.path.split(results_group)
        group = self.data.createGroup(head, tail)
        stations_description = {'s%d' % u: tables.BoolCol() for u in
                                self.stations}

        description = self.reconstruction_description
        description.update(stations_description)
        self.reconstruction = self.data.createTable(group,
            'reconstructions', description)

        description = self.reconstruction_coincidence_description
        description.update(stations_description)
        self.reconstruction_coincidences = \
            self.data.createTable(group, 'coincidences', description)

        return group

    def reconstruct_angles(self, coincidences_group):
        coincidences_group = self.data.getNode(coincidences_group)
        self.data_group = coincidences_group
        coincidences = coincidences_group.coincidences

        self.cluster = coincidences_group._v_attrs.cluster
        self.results_group._v_attrs.cluster = self.cluster

        progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                           pb.ETA()])
        sel_coincidences = coincidences.readWhere('N >= 3')
        for coincidence in progress(sel_coincidences):
            self.reconstruct_individual_stations(coincidence)
            self.reconstruct_cluster_stations(coincidence)

        self.data.flush()

    def reconstruct_individual_stations(self, coincidence):
        coinc_id = coincidence['id']
        event_indexes = self.data_group.c_index[coinc_id]
        events = self.data_group.observables.readCoordinates(event_indexes)

        for event in events:
            if min(event['n1'], event['n3'], event['n4']) >= self.min_n134:
                if self.detector_offsets:
                    theta, phi = self.reconstruct_angle(event, self.detector_offsets[event['station_id']])
                else:
                    theta, phi = self.reconstruct_angle(event)

                if not isnan(theta) and not isnan(phi):
                    self.store_reconstructed_event_from_single_station(coincidence, event, theta, phi)

    def reconstruct_cluster_stations(self, coincidence):
        coinc_id = coincidence['id']
        event_indexes = self.data_group.c_index[coinc_id]
        events = self.data_group.observables.readCoordinates(event_indexes)

        indexes = range(len(events))
        for index_group in itertools.combinations(indexes, 3):
            theta, phi = self.reconstruct_cluster_angle(events, index_group)

            if not isnan(theta) and not isnan(phi):
                self.store_reconstructed_event_from_cluster(coincidence, events, index_group, theta, phi)

    def reconstruct_angle(self, event, offsets=None):
        """Reconstruct angles from a single event"""

        c = 3.00e+8

        if offsets is not None:
            self._correct_offsets(event, offsets)

        dt1 = event['t1'] - event['t3']
        dt2 = event['t1'] - event['t4']

        station = self.cluster.stations[event['station_id']]
        r1, phi1 = station.calc_r_and_phi_for_detectors(1, 3)
        r2, phi2 = station.calc_r_and_phi_for_detectors(1, 4)

        phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                      (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
        theta1 = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
        theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))

        e1 = sqrt(self.rel_theta1_errorsq(theta1, phi, phi1, phi2, r1, r2))
        e2 = sqrt(self.rel_theta2_errorsq(theta2, phi, phi1, phi2, r1, r2))

        theta_wgt = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

        if theta_wgt < 0:
            theta_wgt *= -1
            phi += pi
            phi = (phi + pi) % (2 * pi) - pi

        return theta_wgt, phi

    def reconstruct_cluster_angle(self, events, index_group):
        """Reconstruct angles from a single event"""

        c = 3.00e+8

        t = []
        stations = []
        for event in events:
            timestamp = int(event['ext_timestamp'])
            station = event['station_id']
            offset = int(round(self.station_offsets[station]))
            t.append(timestamp - offset)
            stations.append(station)

        dt1 = t[0] - t[1]
        dt2 = t[0] - t[2]

        r1, phi1 = self.cluster.calc_r_and_phi_for_stations(stations[0], stations[1])
        r2, phi2 = self.cluster.calc_r_and_phi_for_stations(stations[0], stations[2])

        phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                      (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
        theta1 = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
        theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))

        e1 = sqrt(self.rel_theta1_errorsq(theta1, phi, phi1, phi2, r1, r2))
        e2 = sqrt(self.rel_theta2_errorsq(theta2, phi, phi1, phi2, r1, r2))

        theta_wgt = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

        if theta_wgt < 0:
            theta_wgt *= -1
            phi += pi
            phi = (phi + pi) % (2 * pi) - pi

        return theta_wgt, phi

    def store_reconstructed_event_from_single_station(self, coincidence, event,
                                                      reconstructed_theta,
                                                      reconstructed_phi):
        dst_row = self.results_group.reconstructions.row

        dst_row['coinc_id'] = coincidence['id']
        dst_row['N'] = 1
        dst_row['reconstructed_theta'] = reconstructed_theta
        dst_row['reconstructed_phi'] = reconstructed_phi
        dst_row['min_n134'] = min(event['n1'], event['n3'], event['n4'])
        station = self.stations[event['station_id']]
        dst_row['s%d' % station] = True
        dst_row.append()

    def store_reconstructed_event_from_cluster(self, coincidence, events,
                                               index_group, reconstructed_theta,
                                               reconstructed_phi):
        dst_row = self.results_group.reconstructions.row

        dst_row['coinc_id'] = coincidence['id']
        dst_row['N'] = 3
        dst_row['reconstructed_theta'] = reconstructed_theta
        dst_row['reconstructed_phi'] = reconstructed_phi
        for index in index_group:
            station_id = events[index]['station_id']
            station = self.stations[station_id]
            dst_row['s%d' % station] = True
        dst_row.append()


if __name__ == '__main__':
    np.seterr(divide='ignore', invalid='ignore')

    master = Master('master.h5')
    master.main()
