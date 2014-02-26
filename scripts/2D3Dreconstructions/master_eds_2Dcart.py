from math import sqrt
import datetime
import operator
import os.path
import itertools

import tables
import numpy as np
from numpy import arctan2, cos, sin, arcsin, arccos, isnan, pi, linspace
import progressbar as pb
from scipy.optimize import curve_fit
from scipy.stats import norm

from sapphire.esd import download_data
import sapphire.analysis.coincidences
from sapphire.analysis.process_events import (ProcessEvents,
                                              ProcessIndexedEventsWithLINT)
from sapphire.analysis.direction_reconstruction import \
    DirectionReconstruction
from sapphire.analysis.core_reconstruction import (CoreReconstruction,
                                                   CorePositionSolver)
from sapphire.simulations import ldf
from sapphire import storage, clusters


class Master:
    stations = [501, 502, 503, 504, 505, 506, 508]
    datetimerange = (datetime.datetime(2013, 7, 10),
                     datetime.datetime(2013, 7, 11))

    def __init__(self, data_path):
        self.data = tables.openFile(data_path, 'a')

        self.station_groups = ['/s%d' % u for u in self.stations]
        self.cluster = clusters.ScienceParkCluster(self.stations)

        self.detector_offsets = []
        self.station_offsets = []

    def main(self):
        self.download_data()
        self.clean_data()

        self.search_coincidences(window=10000, overwrite=False)

        self.check_coincidences_for_duplicates()

        self.determine_detector_offsets(overwrite=False)
        self.determine_station_offsets(overwrite=False)

        self.reconstruct_direction(overwrite=True)
        self.reconstruct_core_position(overwrite=False)

    def download_data(self):
        start, end = self.datetimerange

        for station, group_path in zip(self.stations, self.station_groups):
            if not group_path in self.data:
                print "Downloading data for station", station
                download_data(self.data, group_path, station,
                              start, end)

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
        self.data.flush()

    def check_coincidences_for_duplicates(self):
        num_stations = len(self.stations)
        coincidences = self.data.root.coincidences

        for station_idx in range(num_stations):
            print "Checking duplicates in coincidences for station %d..." \
                % self.stations[station_idx],
            timestamps = coincidences.observables.readWhere(
                'station_id==station_idx', field='ext_timestamp')
            try:
                assert len(timestamps) == len(set(timestamps))
            except AssertionError:
                raise
            else:
                print 'ok'

        c_index = coincidences.c_index.read()
        print "Checking for duplicate stations in coincidences...",
        found_duplicates = False
        for idx, coords in enumerate(c_index):
            station_ids = \
                coincidences.observables.readCoordinates(
                    coords, field='station_id')
            if len(station_ids) != len(set(station_ids)):
                print "Found duplicate, flushing..."
                found_duplicates = True
                c_index[idx] = []
                coincidences.coincidences.modifyColumn(idx, colname='N', column=0)
        if found_duplicates:
            coincidences.c_index.remove()
            self.data.createVLArray('/coincidences', 'c_index', tables.UInt32Atom())
            for row in c_index:
                coincidences.c_index.append(row)
        else:
            print 'ok'

    def search_coincidences(self, window=200000, overwrite=False):
        if '/coincidences' not in self.data or overwrite:
            print "Searching for coincidences..."
            coincidences = sapphire.analysis.coincidences.Coincidences(
                self.data, '/coincidences', self.station_groups,
                overwrite=overwrite)
            coincidences.search_coincidences(window=window) # Add this: window=10000
            coincidences.process_events()
            coincidences.store_coincidences(self.cluster)

    def reconstruct_direction(self, overwrite=False):
        if not '/reconstructions' in self.data or overwrite:
            print "Reconstructing direction..."
            reconstruction = ClusterDirectionReconstruction(self.data,
                                self.stations, '/reconstructions',
                                detector_offsets=self.detector_offsets,
                                station_offsets=self.station_offsets,
                                overwrite=overwrite)
            reconstruction.reconstruct_angles('/coincidences')

    def reconstruct_core_position(self, overwrite=False):
        if not '/core_reconstructions' in self.data or overwrite:
            print "Reconstructing core position..."
            reconstruction = CoreReconstruction(self.data, self.stations,
                                '/core_reconstructions',
                                overwrite=overwrite)
            reconstruction.reconstruct_core_positions('/coincidences')

    def determine_detector_offsets(self, overwrite=False):
        offsets_group = '/detector_offsets'
        if offsets_group in self.data and not overwrite:
            self.detector_offsets = self.data.root.detector_offsets.read()
        else:
            for station_id, station_group in enumerate(self.station_groups):
                process = ProcessEvents(self.data, station_group)
                offsets = process.determine_detector_timing_offsets()
                print "Offsets for station %d: %s" % (station_id, offsets)
                self.detector_offsets.append(offsets)
            if offsets_group in self.data:
                self.data.removeNode(offsets_group)
            group, node = os.path.split(offsets_group)
            self.data.createArray(group, node, self.detector_offsets)

    def determine_station_offsets(self, overwrite=False):
        offsets_group = '/station_offsets'
        if offsets_group in self.data and not overwrite:
            self.station_offsets = self.data.root.station_offsets.read()
        else:
            ref_group = '/s501'
            station_groups = list(self.station_groups)
            station_groups.remove(ref_group)

            gauss = lambda x, N, mu, sigma: N * norm.pdf(x, mu, sigma)
            bins = linspace(-1e3, 1e3, 101)

            for station_id, station_group in enumerate(station_groups):
                coincidences = sapphire.analysis.coincidences.Coincidences(
                    self.data, coincidence_group=None,
                    station_groups=[ref_group, station_group])
                c_index, timestamps = coincidences._search_coincidences()

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

            if offsets_group in self.data:
                self.data.removeNode(offsets_group)
            group, node = os.path.split(offsets_group)
            self.data.createArray(group, node, self.station_offsets)


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
            sel_events = [events[u] for u in index_group]
            theta, phi = self.reconstruct_cluster_angle(sel_events)

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

    def reconstruct_cluster_angle(self, events):
        """Reconstruct angles from three events"""

        c = 3.00e+8

        t = []
        stations = []
        for event in events:
            timestamp = int(event['ext_timestamp'])
            station = event['station_id']
            arrival_times = [event[u] for u in 't1', 't2', 't3', 't4' if
                             event[u] != -999.]
            arrival_times.sort()
            # FIXME: should check for three low condition (< 1% ?)
            if len(arrival_times) >= 2:
                trigger_offset = arrival_times[1] - arrival_times[0]
            else:
                trigger_offset = 0.
            offset = self.station_offsets[station]
            # FIXME: ext_timestamp (long) + float loses precision
            try:
                correction = int(round(offset + trigger_offset))
            except:
                print arrival_times
                correction = int(round(offset))
            t.append(timestamp - correction)
            stations.append(station)

        dt1 = t[1] - t[0]
        dt2 = t[2] - t[0]
        dt1 = dt1 * 1e-9
        dt2 = dt2 * 1e-9

        x0, y0, z0, _ = self.cluster.stations[stations[0]].get_xyzalpha_coordinates()
        x1, y1, z1, _ = self.cluster.stations[stations[1]].get_xyzalpha_coordinates()
        x2, y2, z2, _ = self.cluster.stations[stations[2]].get_xyzalpha_coordinates()

        dx1 = x1 - x0
        dy1 = y1 - y0

        dx2 = x2 - x0
        dy2 = y2 - y0

        ux = c * dt2 * dx1 - c * dt1 * dx2
        uy = c * dt2 * dy1 - c * dt1 * dy2

        vz = dx1 * dy2 - dx2 * dy1

        usquared = ux * ux + uy * uy
        vzsquared = vz * vz


        phi = arctan2(-ux * vz,uy * vz)
        theta = arcsin(sqrt(usquared/vzsquared))

        if isnan(theta):
            phi = np.nan

        return theta, phi



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

    master = Master('master_2Dcart_1dag.h5')
    master.main()
