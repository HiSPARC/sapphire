""" Search for coincidences between HiSPARC stations

    This module can be used to search for coincidences between several
    HiSPARC stations. To skip this and directly download coincidences
    use :func:`~sapphire.esd.download_coincidences`, this is slightly
    less flexible because you can not choose the coincidence window.

    For regular usage, download events from the ESD and use the
    :class:`CoincidencesESD` class. Example usage::

        import datetime

        import tables

        from sapphire import CoincidencesESD, download_data

        STATIONS = [501, 503, 506]
        START = datetime.datetime(2013, 1, 1)
        END = datetime.datetime(2013, 1, 2)


        if __name__ == '__main__':
            station_groups = ['/s%d' % u for u in STATIONS]

            data = tables.open_file('data.h5', 'w')
            for station, group in zip(STATIONS, station_groups):
                download_data(data, group, station, START, END)

            coin = CoincidencesESD(data, '/coincidences', station_groups)
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

    .. note::
        For better compatibility with other modules, such as
        :mod:`~sapphire.analysis.reconstructions`, it is recommended to use the
        subclass :class:`CoincidencesESD` instead. This is an old class with a
        different way to store the results. It can be used, however, on *raw*
        data downloaded from the public database.

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

    Once the coincidences are stored, there will be a `coincidences` table in
    the group. This table has multiple columns used for storing simulation
    inputs like shower direction and energy. At this point, they contain no
    information. They are only useful if the original event tables were created
    by simulations, instead of real detector data. The useful columns are:

        * ``id``: the index of this coincidence. This index is identical for
          the ``coincidence`` and ``c_index`` tables.
        * ``timestamp``: the timestamp of the event in seconds
        * ``nanoseconds``: the subsecond part of the timestamp in nanoseconds
        * ``ext_timestamp``: the timestamp of the event in nanoseconds (equal
          to timestamp * 1000000000 + nanoseconds)
        * ``N``: the number of stations participating in this coincidence
        * ``s0``, ``s1``, ...: whether the first (0), second (1) or other
          stations participated in the coincidence.

    The coincidences group furthermore contains the table ``c_index`` to track
    down the individual events making up the coincidence. The ``c_index`` table
    gives indexes for the individual events inside the original event tables.

    If you have obtained a particular coincidence from the ``coincidences``
    table, the ``id`` is the index into all these tables. For example, looking
    up the source events making up the 40th coincidence::

        >>> group = data.root.coincidences
        >>> idx = 40
        >>> group.coincidences[idx]

    can be done in the following way (each row in the ``c_index`` table is a
    (station index, event index) pair)::

        >>> for event_idx in group.c_index[idx]:
        ...     event = group.observables[event_idx]

    The ``event`` is one of the source events, processed to determine particle
    arrival times from the raw traces (if available). It has a
    ``station_id`` attribute which is an index into ``station_groups``.

    """

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

    def search_and_store_coincidences(self, window=10000):
        """Search, process and store coincidences.

        This is a semi-automatic method to search for coincidences,
        process the events making up the coincidences and then store the
        results in the coincidences group.

        If you want to make use of non-default parameters like time
        shifts or overwriting previously processed events, please call
        the individual methods.  See the class docstring.

        """
        self.search_coincidences(window=window)
        self.process_events()
        self.store_coincidences()

    def search_coincidences(self, window=10000, shifts=None, limit=None):
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

        :param window: the coincidence time window in nanoseconds. All events
            with delta t's smaller than this window will be considered a
            coincidence.
        :param shifts: optionally shift a station's data in time.  This
            can be useful if a station has a misconfigured GPS clock.
            Expects a list of shifts, one for each station, in seconds.
            Use 'None' for no shift.
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

        if len(c_index) == 0:
            return

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
                if self.progress:
                    print "Processing coincidence events with traces"
                Process = process_events.ProcessIndexedEventsWithLINT
            else:
                if self.progress:
                    print "Processing coincidence events without traces"
                Process = process_events.ProcessIndexedEventsWithoutTraces

            process = Process(self.data, station_group, index,
                              progress=self.progress)
            process.process_and_store_results(overwrite=overwrite)

    def store_coincidences(self):
        """Store the previously found coincidences.

        After you have searched for coincidences, you can store the
        more user-friendly results in the coincidences group using this
        method.

        """
        self.c_index = []
        self.coincidences = self.data.create_table(self.coincidence_group,
                                                   'coincidences',
                                                   storage.Coincidence)
        self.observables = self.data.create_table(self.coincidence_group,
                                                  'observables',
                                                  storage.EventObservables)

        for coincidence in pbar(self.coincidence_group._src_c_index,
                                show=self.progress):
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
            event_index = int(event_desc[2])

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

    def _search_coincidences(self, window=10000, shifts=None, limit=None):
        """Search for coincidences

        Search for coincidences in a set of PyTables event tables, optionally
        shifting the data in time.  This is necessary when one wants to
        compare the timestamps of stations who use a different time (as in
        GPS, UTC or local time).  This function searches for events which
        occured almost at the same time and thus might be the result of an
        extended air shower.

        :param window: the time window in nanoseconds which will be searched
            for coincidences.  Events falling outside this window will not be
            part of the coincidence.  Default: 10000 (i.e. 10 us).
        :param shifts: a list of time shifts in seconds, use 'None' for no
            shift.
        :param limit: limit the number of events which are processed.

        :return: coincidences, timestamps. First a list of coincidences, which
            each consist of a list with indexes into the timestamps array as a
            pointer to the events making up the coincidence. Then, a list of
            tuples.  Each tuple consists of a timestamp followed by an index
            into the stations list which designates the detector
            station which measured the event, and finally an index into that
            station's event table.

        """
        # get the 'events' tables from the groups or groupnames
        event_tables = []
        for station_group in self.station_groups:
            station_group = self.data.get_node(station_group)
            if 'events' in station_group:
                event_tables.append(self.data.get_node(station_group,
                                                       'events'))

        timestamps = self._retrieve_timestamps(event_tables, shifts, limit)
        coincidences = self._do_search_coincidences(timestamps, window)

        return coincidences, timestamps

    def _retrieve_timestamps(self, event_tables, shifts=None, limit=None):
        """Retrieve all timestamps from all stations, optionally shifting them

        This function retrieves the timestamps from a list of event tables and
        will optionally shift the timestamps by a given amount.  This is
        necessary when one wants to compare the timestamps of stations who use
        a different time (as in GPS, UTC or local time).

        :param event_tables: a list of HiSPARC event tables, usually from
            different stations.
        :param shifts: a list of time shifts in seconds, use 'None' for no
            shift.
        :param limit: limit the number of events which are processed.

        :return: list of tuples.  Each tuple consists of a timestamp followed
            by an index into the stations list which designates the detector
            station which measured the event, and finally an index of the
            event into the station's event table.

        """
        # calculate the shifts in nanoseconds and cast them to int.
        # (prevent upcasting timestamps to float64 further on)
        if shifts is not None:
            shifts = [int(shift * 1e9) if shift is not None else shift
                      for shift in shifts]

        timestamps = []
        for s_id, event_table in enumerate(event_tables):
            ts = [(x, s_id, j) for j, x in
                  enumerate(event_table.col('ext_timestamp')[:limit])]
            try:
                # shift data. carefully avoid upcasting (we're adding two
                # ints, which is an int, and casting that back to uint64. if
                # we're not careful, an intermediate value will be a float64,
                # which doesn't hold the precision to store nanoseconds.
                ts = [(np.uint64(int(x) + shifts[i]), i, j)
                      for x, i, j in ts]
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
    """Store coincidences specifically using the ESD

    This is a subclass of :class:`Coincidences`. This subclass stores the paths
    to the station_groups that where used to look for coincidences in a
    lookup-table. The c_index stores the station and event id for each event in
    the coincidence, allowing you to find the original event row.

    Suppose you want to search for coincidences between stations 501 and 503.
    First, download the data for these stations from the ESD. Suppose you
    stored the data in the '/s501' and '/s503' groups in the file 'data'.
    Then::

        >>> station_groups = ['/s501', '/s503']
        >>> coin = CoincidencesESD(data, '/coincidences', station_groups)
        >>> coin.search_and_store_coincidences()

    If you want a more manual method, replace the last line with::

        >>> coin.search_coincidences(window=5000)
        >>> coin.process_events()
        >>> coin.store_coincidences()

    You can then provide different parameters to the individual methods.
    See the corresponding docstrings.

    Once the coincidences are stored, there will be a `coincidences` table in
    the group. This table has multiple columns used for storing simulation
    inputs like shower direction and energy. At this point, they contain no
    information. They are only used if the event and coincidence tables were
    created by simulations, instead of real detector data. The useful columns
    are:

        * ``id``: the index of this coincidence. This index is identical for
          the ``coincidence`` and ``c_index`` tables.
        * ``timestamp``: the timestamp of the event in seconds
        * ``nanoseconds``: the subsecond part of the timestamp in nanoseconds
        * ``ext_timestamp``: the timestamp of the event in nanoseconds (equal
          to timestamp * 1000000000 + nanoseconds)
        * ``N``: the number of stations participating in this coincidence
        * ``s0``, ``s1``, ...: for each station indicate whether it
          participated in the coincidence.

    The coincidences group furthermore contains the tables ``s_index`` and
    ``c_index`` to track down the individual events making up the coincidence.
    The ``s_index`` table contains a row for each station, pointing to the
    station's event tables. The ``c_index`` table gives indexes for the
    individual events inside those tables.

    If you have obtained a particular coincidence from the ``coincidences``
    table, the ``id`` is the index into all these tables. For example, looking
    up the source events making up the 40th coincidence::

        >>> group = data.root.coincidences
        >>> idx = 40
        >>> group.coincidences[idx]

    can be done in the following way (each row in the ``c_index`` table is a
    (station index, event index) pair)::

        >>> for station_idx, event_idx in group.c_index[idx]:
        ...     station_path = group.s_index[station_idx]
        ...     event_group = data.get_node(station_path, 'events')
        ...     event = event_group[event_idx]

    The ``event`` is one of the events in the coincidence.

    """
    def search_and_store_coincidences(self, window=10000,
                                      station_numbers=None):
        """Search and store coincidences.

        This is a semi-automatic method to search for coincidences
        and then store the results in the coincidences group.

        """
        self.search_coincidences(window=window)
        self.store_coincidences(station_numbers=station_numbers)

    def search_coincidences(self, window=10000, shifts=None, limit=None):
        """Search for coincidences.

        Search all data in the station_groups for coincidences, and store
        rudimentary coincidence data in attributes.  This data might be useful,
        but is very basic.  You can call the :meth:`store_coincidences` method
        to store the coincidences in an easier format in the coincidences
        group.

        If you want to process the preliminary results: they are stored in the
        attributes ``_src_c_index`` and ``_src_timestamps``.  The
        former is a list of coincidences, which each consist of a list with
        indexes into the timestamps array as a pointer to the events making up
        the coincidence. The latter is a list of tuples.  Each tuple consists
        of a timestamp followed by an index into the stations list which
        designates the detector station which measured the event, and finally
        an index into that station's event table.

        :param window: the coincidence time window.  All events with delta
            t's smaller than this window will be considered a coincidence.
        :param shifts: optionally shift a station's data in time.  This
            can be useful if a station has a misconfigured GPS clock.
            Expects a list of shifts, one for each station.
        :param limit: optionally limit the search for this number of
            events.

        """
        c_index, timestamps = self._search_coincidences(window, shifts, limit)
        self._src_timestamps = timestamps
        self._src_c_index = c_index

    def store_coincidences(self, station_numbers=None):
        """Store the previously found coincidences.

        After having searched for coincidences, you can store the more
        user-friendly results in the ``coincidences`` group using this
        method. It also created a ``c_index`` and ``s_index`` table to
        find the source events.

        :param station_numbers: optional list of station_numbers.
            If given these will be used to attach correct numbers to the
            station column names in the coincidences table. Otherwise
            they will simply be numbered by id. This list must be the
            same length as the station_groups.

        """
        n_coincidences = len(self._src_c_index)
        if station_numbers is not None:
            if len(station_numbers) != len(self.station_groups):
                raise RuntimeError(
                    "Number of station numbers must equal number of groups.")
            self.station_numbers = station_numbers
            s_columns = {'s%d' % number: tables.BoolCol(pos=p)
                         for p, number in enumerate(station_numbers, 12)}
        else:
            self.station_numbers = None
            s_columns = {'s%d' % n: tables.BoolCol(pos=(n + 12))
                         for n, _ in enumerate(self.station_groups)}

        description = storage.Coincidence
        description.columns.update(s_columns)
        self.coincidences = self.data.create_table(
            self.coincidence_group, 'coincidences', description,
            expectedrows=n_coincidences)

        self.c_index = []

        for coincidence in pbar(self._src_c_index, show=self.progress):
            self._store_coincidence(coincidence)

        c_index = self.data.create_vlarray(
            self.coincidence_group, 'c_index', tables.UInt32Col(shape=2),
            expectedrows=n_coincidences)
        for observables_idx in pbar(self.c_index, show=self.progress):
            c_index.append(observables_idx)
        c_index.flush()

        s_index = self.data.create_vlarray(
            self.coincidence_group, 's_index', tables.VLStringAtom(),
            expectedrows=len(self.station_groups))
        for station_group in self.station_groups:
            s_index.append(station_group)
        s_index.flush()

    def _store_coincidence(self, coincidence):
        """Store a single coincidence in the coincidence group.

        Stores coincidence in the coincidences table and references
        to the observables making up each coincidence in ``c_index``.

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
            if self.station_numbers is not None:
                station_number = self.station_numbers[station_id]
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
        :meth:`~Coincidences.search_coincidences`.
    :param timestamps: the list of timestamps, as returned by
        :meth:`~Coincidences.search_coincidences`.
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
