""" Search for coincidences between HiSPARC stations

    This module can be used to search for coincidences between several
    HiSPARC stations. To skip this and directly download coincidences
    use :func:`~sapphire.esd.download_coincidences`, this is slightly
    less flexible because you can not choose the coincidence window.

    Example usage::

        import datetime

        import tables

        from sapphire.publicdb import download_data
        from sapphire.analysis import coincidences


        STATIONS = [501, 503, 506]
        START = datetime.datetime(2013, 1, 1)
        END = datetime.datetime(2013, 1, 2)


        if __name__ == '__main__':
            station_groups = ['/s%d' % u for u in STATIONS]

            data = tables.open_file('data.h5', 'w')
            for station, group in zip(STATIONS, station_groups):
                download_data(data, group, station, START, END)

            coin = coincidences.CoincidencesESD(data, '/coincidences',
                                                station_groups)
            coin.search_and_store_coincidences()

"""
import os.path

import tables
import numpy as np
from progressbar import ProgressBar, ETA, Bar, Percentage

from . import process_events
from .. import storage
from ..utils import pbar


class Coincidences(object):
    """Search for and store coincidences between HiSPARC stations.

    Suppose you want to search for coincidences between stations 501 and
    503.  First, download the data for these stations (with, or without
    traces, depending on your intentions).  Suppose you stored the data in
    the '/s501' and '/s503' groups in the 'data' file.  Then::

        >>> station_groups = ['/s501', '/s503']
        >>> coin = Coincidences(data, '/coincidences', station_groups)
        >>> coin.search_and_store_coincidences()

    If you want a more manual method, replace the last line with::

        >>> coin.search_coincidences(window=50000)
        >>> coin.process_events()
        >>> coin.store_coincidences()

    You can then provide different parameters to the individual methods.
    See the corresponding docstrings.

    .. note::
        For better compatibility with other modules, such as
        :mod:`~sapphire.analysis.reconstructions`, it is recommended
        to use the subclass :class:`CoincidencesESD` instead.

    """
    ProcessWithTraces = process_events.ProcessIndexedEventsWithLINT
    ProcessWithoutTraces = process_events.ProcessIndexedEventsWithoutTraces

    def __init__(self, data, coincidence_group, station_groups,
                 overwrite=False, progress=True):
        """Initialize the class.

        :param data: the PyTables datafile.
        :param coincidence_group: the destination group.
        :param station_groups: a list of groups containing the station
            data.
        :param overwrite: if True, overwrite a previous coincidences
            group.
        :param progress: if True, show a progressbar while storing
            coincidences.

        """
        self.data = data
        if coincidence_group is not None:
            if coincidence_group in self.data:
                if overwrite:
                    self.data.remove_node(coincidence_group, recursive=True)
                else:
                    raise RuntimeError("Group %s already exists in datafile, "
                                       "and overwrite is False" %
                                       coincidence_group)
            head, tail = os.path.split(coincidence_group)
            self.coincidence_group = data.create_group(head, tail,
                                                       createparents=True)
        self.station_groups = station_groups

        self.trig_threshold = .5
        self.overwrite = overwrite
        self.progress = progress

    def search_and_store_coincidences(self, cluster=None):
        """Search, process and store coincidences.

        This is a semi-automatic method to search for coincidences,
        process the events making up the coincidences and then store the
        results in the coincidences group.

        If you want to make use of non-default parameters like coincidence
        window lenghts, time shifts or overwriting previously processed
        events, please call the individual methods.  See the class
        docstring.

        """
        self.search_coincidences()
        self.process_events()
        self.store_coincidences(cluster)

    def search_coincidences(self, window=200000, shifts=None, limit=None):
        """Search for coincidences.

        Search all data in the station_groups for coincidences, and store
        rudimentary coincidence data in the coincidences group.  This data
        might be useful, but is very basic.  You can call the
        :meth:`store_coincidences` method to store the coincidences in an
        easier format in the coincidences group.

        If you want to process the preliminary results: they are stored in
        _src_c_index and _src_timestamps.  The former is a list of
        coincidences, which each consist of a list with indexes into the
        timestamps array as a pointer to the events making up the
        coincidence. The latter is a list of tuples.  Each tuple consists
        of a timestamp followed by an index into the stations list which
        designates the detector station which measured the event, and
        finally an index into that station's event table.

        :param window: the coincidence time window.  All events with delta
            t's smaller than this window will be considered a coincidence.
        :param shifts: optionally shift a station's data in time.  This
            can be useful if a station has a misconfigured GPS clock.
            Expects a list of shifts for each station.
        :param limit: optionally limit the search for this number of
            events.

        """
        c_index, timestamps = \
            self._search_coincidences(window, shifts, limit)
        timestamps = np.array(timestamps, dtype=np.uint64)
        self.data.create_array(self.coincidence_group, '_src_timestamps',
                               timestamps)
        src_c_index = self.data.create_vlarray(self.coincidence_group,
                                               '_src_c_index',
                                               tables.UInt32Atom())
        for coincidence in c_index:
            src_c_index.append(coincidence)

    def process_events(self, overwrite=None):
        """Process events using :mod:`~sapphire.analysis.process_events`

        Events making up the coincidences are processed to obtain
        observables like number of particles and particle arrival times.

        :param overwrite: if True, overwrite the events tables in the
            station groups.

        """
        if overwrite is None:
            overwrite = self.overwrite

        c_index = self.coincidence_group._src_c_index.read()
        timestamps = self.coincidence_group._src_timestamps.read()

        selected_timestamps = []
        for coincidence in c_index:
            for event in coincidence:
                selected_timestamps.append(timestamps[event])
        full_index = np.array(selected_timestamps)

        for station_id, station_group in enumerate(self.station_groups):
            station_group = self.data.get_node(station_group)
            selected = full_index.compress(full_index[:, 1] == station_id,
                                           axis=0)
            index = selected[:, 2]

            if 'blobs' in station_group:
                print "Processing coincidence events with traces"
                Process = self.ProcessWithTraces
            else:
                print "Processing coincidence events without traces"
                Process = self.ProcessWithoutTraces

            process = Process(self.data, station_group, index)
            process.process_and_store_results(overwrite=overwrite)

    def store_coincidences(self, cluster=None):
        """Store the previously found coincidences.

        After you have searched for coincidences, you can store the
        more user-friendly results in the coincidences group using this
        method.

        :param cluster: optionally store a
            :class:`~sapphire.clusters.BaseCluster` instance in the
            coincidences group for future reference.

        """
        if cluster:
            self.coincidence_group._v_attrs.cluster = cluster

        self.c_index = []
        self.coincidences = self.data.create_table(self.coincidence_group,
                                                   'coincidences',
                                                   storage.Coincidence)
        self.observables = self.data.create_table(self.coincidence_group,
                                                  'observables',
                                                  storage.EventObservables)

        # ProgressBar does not work for empty iterables.
        if len(self.coincidence_group._src_c_index):
            src_c_index = pbar(self.coincidence_group._src_c_index,
                               show=self.progress)
            for coincidence in src_c_index:
                self._store_coincidence(coincidence)
        else:
            print "Creating empty tables, no coincidences found"
            for coincidence in self.coincidence_group._src_c_index:
                self._store_coincidence(coincidence)

        c_index = self.data.create_vlarray(self.coincidence_group, 'c_index',
                                           tables.UInt32Col())
        for coincidence in self.c_index:
            c_index.append(coincidence)
        c_index.flush()
        self.c_index = c_index

    def _store_coincidence(self, coincidence):
        """Store a single coincidence in the coincidence group.

        Stores in the coincidences table, the c_index, and the individual
        events in the observables table.

        """
        row = self.coincidences.row
        coincidence_id = len(self.coincidences)
        row['id'] = coincidence_id
        row['N'] = len(coincidence)

        observables_idx = []
        timestamps = []
        for index in coincidence:
            event_desc = self.coincidence_group._src_timestamps[index]
            station_id = event_desc[1]
            event_index = event_desc[2]

            group = self.data.get_node(self.station_groups[station_id])
            event = group.events[event_index]
            idx = self._store_event_in_observables(event, coincidence_id,
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

    def _store_event_in_observables(self, event, coincidence_id,
                                    station_id):
        """Store a single event in the observables table."""

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

    def _search_coincidences(self, window=200000, shifts=None, limit=None):
        """Search for coincidences

        Search for coincidences in a set of PyTables event tables, optionally
        shifting the data in time.  This is necessary when one wants to
        compare the timestamps of stations who use a different time (as in
        GPS, UTC or local time).  This function searches for events which
        occured almost at the same time and thus might be the result of an
        extended air shower.

        :param data: the PyTables data file
        :param stations: a list of HiSPARC event tables (normally from
            different stations, hence the name)
        :param window: the time window in nanoseconds which will be searched
            for coincidences.  Events falling outside this window will not be
            part of the coincidence.  Default: 200000 (i.e. 200 us).
        :param shifts: a list of time shifts which may contain 'None'.
        :param limit: limit the number of events which are processed

        :return: coincidences, timestamps. First a list of coincidences, which
            each consist of a list with indexes into the timestamps array as a
            pointer to the events making up the coincidence. Then, a list of
            tuples.  Each tuple consists of a timestamp followed by an index
            into the stations list which designates the detector
            station which measured the event, and finally an index into that
            station's event table.

        Example usage::

            >>> import tables
            >>> from hisparc.analysis.coincidences import search_coincidences
            >>> data = tables.open_file('test.h5', 'a')
            >>> coincidences, timestamps = search_coincidences(data,
            ... ['/hisparc/station501', '/hisparc/station502',
            ... '/hisparc/station503', '/hisparc/station504',
            ... '/hisparc/station505'], shifts=[None, None, -15, None, None])
            >>> coincidences[:3]
            [[73, 74], [81, 82], [98, 99]]
            >>> timestamps[73], timestamps[74]
            ((1235433610410730837, 0, 23), (1235433610410731004, 2, 17))

        """
        # get the 'events' tables from the groups or groupnames
        event_tables = []
        for station_group in self.station_groups:
            station_group = self.data.get_node(station_group)
            if 'events' in station_group:
                event_tables.append(self.data.get_node(station_group,
                                                       'events'))
        stations = event_tables

        # calculate the shifts in nanoseconds and cast them to long.
        # (prevent upcasting timestamps to float64 further on)
        if shifts:
            for i in range(len(shifts)):
                if shifts[i]:
                    shifts[i] = int(shifts[i] * 1e9)

        timestamps = self._retrieve_timestamps(stations, shifts, limit)
        coincidences = self._do_search_coincidences(timestamps, window)

        return coincidences, timestamps

    def _retrieve_timestamps(self, stations, shifts=None, limit=None):
        """Retrieve all timestamps from all stations, optionally shifting them

        This function retrieves the timestamps from a list of event tables and
        will optionally shift the timestamps by a given amount.  This is
        necessary when one wants to compare the timestamps of stations who use
        a different time (as in GPS, UTC or local time).

        :param stations: a list of HiSPARC event tables (normally from
            different stations, hence the name)
        :param shifts: a list of time shifts which may contain 'None'.
        :param limit: limit the number of events which are processed

        :return: list of tuples.  Each tuple consists of a timestamp followed
            by an index into the stations list which designates the detector
            station which measured the event, and finally an index of the
            event into the station's event table.

        """
        timestamps = []
        for i in range(len(stations)):
            ts = [(x['ext_timestamp'], i, j) for j, x in
                  enumerate(stations[i][:limit])]
            try:
                # shift data. carefully avoid upcasting (we're adding two
                # longs, which is a long, and casting that back to uint64. if
                # we're not careful, an intermediate value will be a float64,
                # which doesn't hold the precision to store nanoseconds.
                ts = [(np.uint64(int(x[0]) + shifts[i]), x[1], x[2]) for x in
                      ts]
            except (TypeError, IndexError):
                # shift is None or doesn't exist
                pass
            timestamps.extend(ts)

        # sort the timestamps
        timestamps.sort()

        return timestamps

    def _do_search_coincidences(self, timestamps, window):
        """Search for coincidences in a set of timestamps

        Given a set of timestamps, search for coincidences.  That is, search
        for events which occured almost at the same time and thus might be the
        result of an extended air shower.

        :param timestamps: a list of tuples (timestamp, station_idx,
            event_idx) which will be searched
        :param window: the time window in nanoseconds which will be searched
            for coincidences.  Events falling outside this window will not be
            part of the coincidence.

        :return: a list of coincidences, which each consist of a list with
            indexes into the timestamps array as a pointer to the events
            making up the coincidence

        """
        coincidences = []

        # traverse all timestamps
        prev_coincidence = []

        if self.progress and len(timestamps):
            pbar = ProgressBar(maxval=len(timestamps),
                               widgets=[Percentage(), Bar(), ETA()]).start()

        for i in xrange(len(timestamps)):

            # build coincidence, starting with the current timestamp
            c = [i]
            t0 = timestamps[i][0]

            # traverse the rest of the timestamps
            for j in xrange(i + 1, len(timestamps)):
                # if a timestamp is within the coincidence window, add it
                if timestamps[j][0] - t0 < window:
                    c.append(j)
                else:
                    # coincidence window has passed, break for-loop
                    break

            # if we have more than one event in the coincidence, save it
            if len(c) > 1:
                # is this coincidence part of the previous coincidence?
                is_part_of_prev = np.array([u in prev_coincidence
                                            for u in c]).all()
                if not is_part_of_prev:
                    # no, so it's a new one
                    coincidences.append(c)
                    prev_coincidence = c

            if self.progress and not i % 5000:
                pbar.update(i)

        if self.progress and len(timestamps):
            pbar.finish()

        return coincidences


class CoincidencesESD(Coincidences):
    """Store coincidences differently for the ESD

    This subclass stores the paths to the station_groups that where
    used to look for coincidences in a lookup-table, and also
    stores the original station info and event_id for each coincidence.

    """
    def search_and_store_coincidences(self, cluster=None):
        """Search and store coincidences.

        This is a semi-automatic method to search for coincidences
        and then store the results in the coincidences group.

        """
        self.search_coincidences()
        self.store_coincidences(cluster)

    def search_coincidences(self, window=200000, shifts=None, limit=None):
        """Search for coincidences.

        Instead of storing the results in the tables `_src_c_index` and
        `_src_timestamps`, they are stored in attributes by the same
        name in the class.

        """
        c_index, timestamps = self._search_coincidences(window, shifts, limit)
        self._src_timestamps = np.array(timestamps, dtype=np.uint64)
        self._src_c_index = c_index

    def store_coincidences(self, cluster=None):
        """Store the previously found coincidences.

        After having searched for coincidences, you can store the more
        user-friendly results in the `coincidences` group using this
        method. It also created a `c_index` and `s_index` table to find
        the source events.

        """
        if cluster:
            self.cluster = cluster
            self.coincidence_group._v_attrs.cluster = cluster
            s_columns = {'s%d' % station.number: tables.BoolCol(pos=p)
                         for p, station in enumerate(cluster.stations, 12)}
        else:
            self.cluster = None
            s_columns = {'s%d' % n: tables.BoolCol(pos=(n + 12))
                         for n, _ in enumerate(self.station_groups)}

        self.c_index = []

        description = storage.Coincidence
        description.columns.update(s_columns)
        self.coincidences = self.data.create_table(self.coincidence_group,
                                                   'coincidences', description)

        # ProgressBar does not work for empty iterables.
        if len(self._src_c_index):
            src_c_index = pbar(self._src_c_index, show=self.progress)
            for coincidence in src_c_index:
                self._store_coincidence(coincidence)
        else:
            print "Creating empty tables, no coincidences found"
            for coincidence in self._src_c_index:
                self._store_coincidence(coincidence)

        c_index = self.data.create_vlarray(self.coincidence_group, 'c_index',
                                           tables.UInt32Col(shape=2))
        for coincidence in self.c_index:
            c_index.append(coincidence)
        c_index.flush()
        self.c_index = c_index

        s_index = self.data.create_vlarray(self.coincidence_group, 's_index',
                                           tables.VLStringAtom())
        for station_group in self.station_groups:
            s_index.append(station_group)
        s_index.flush()

    def _store_coincidence(self, coincidence):
        """Store a single coincidence in the coincidence group.

        Stores coincidence in the coincidences table and references
        to the observables making up each coincidence in `c_index`.

        """
        row = self.coincidences.row
        coincidence_id = len(self.coincidences)
        row['id'] = coincidence_id
        row['N'] = len(coincidence)

        observables_idx = []
        timestamps = []
        for index in coincidence:
            event_desc = self._src_timestamps[index]
            station_id = event_desc[1]
            event_index = event_desc[2]
            if self.cluster:
                station_number = self.cluster.stations[station_id].number
                row['s%d' % station_number] = True
            else:
                row['s%d' % station_id] = True

            group = self.data.get_node(self.station_groups[station_id])
            event = group.events[event_index]
            observables_idx.append((station_id, event_index))
            timestamps.append((event['ext_timestamp'], event['timestamp'],
                               event['nanoseconds']))

        first_timestamp = sorted(timestamps)[0]
        row['ext_timestamp'], row['timestamp'], row['nanoseconds'] = \
            first_timestamp
        row.append()
        self.c_index.append(observables_idx)
        self.coincidences.flush()


def get_events(data, stations, coincidence, timestamps, get_raw_traces=False):
    """Get event data of a coincidence

    Return a list of events making up a coincidence.

    :param data: the PyTables data file
    :param stations: a list of HiSPARC event tables (normally from
        different stations, hence the name)
    :param coincidence: a coincidence, as returned by
        :func:`search_coincidences`.
    :param timestamps: the list of timestamps, as returned by
        :func:`search_coincidences`.
    :param get_raw_traces: boolean.  If true, return the compressed ADC
        values instead of the uncompressed traces.

    :return: a list of tuples.  Each tuple consists of (station, event,
        traces), where event is the event row from PyTables and traces is
        a list of the uncompressed traces.

    """
    events = []
    for event in coincidence:
        timestamp, station, index = timestamps[event]
        process = process_events.ProcessEvents(data, stations[station])
        event = process.source[index]
        if not get_raw_traces:
            baseline = np.where(event['baseline'] != -999, event['baseline'],
                                200)[np.where(event['traces'] >= 0)]
            # transpose to get expected format
            traces = (process.get_traces_for_event(event) - baseline).T
            traces = traces * -0.57
        else:
            traces = [process.group.blobs[x] for x in event['traces']]
        events.append((stations[station], event, traces))

    return events
