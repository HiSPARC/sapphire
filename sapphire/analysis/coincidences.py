""" Search for coincidences between HiSPARC stations

    This module can be used to search for coincidences between several
    HiSPARC stations.

    To search for coincidences, use the :func:`search_coincidences`
    function.
"""

import numpy as np
import time


def search_coincidences(data, stations, window=200000, shifts=None, limit=None):
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
        >>> data = tables.openFile('test.h5', 'a')
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
    for station_group in stations:
        station_group = data.getNode(station_group)
        if 'events' in station_group:
            event_tables.append(data.getNode(station_group, 'events'))
    stations = event_tables

    # calculate the shifts in nanoseconds and cast them to long.
    # (prevent upcasting timestamps to float64 further on)
    if shifts:
        for i in range(len(shifts)):
            if shifts[i]:
                shifts[i] = long(shifts[i] * 1e9)

    timestamps = retrieve_timestamps(stations, shifts, limit)
    coincidences = do_search_coincidences(timestamps, window)

    return coincidences, timestamps


def retrieve_timestamps(stations, shifts=None, limit=None):
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
            ts = [(np.uint64(long(x[0]) + shifts[i]), x[1], x[2]) for x in
                  ts]
        except (TypeError, IndexError):
            # shift is None or doesn't exist
            pass
        timestamps.extend(ts)

    # sort the timestamps
    timestamps.sort()

    return timestamps


def do_search_coincidences(timestamps, window=200000):
    """Search for coincidences in a set of timestamps

    Given a set of timestamps, search for coincidences.  That is, search
    for events which occured almost at the same time and thus might be the
    result of an extended air shower.

    :param timestamps: a list of tuples (timestamp, station_idx,
        event_idx) which will be searched
    :param window: the time window in nanoseconds which will be searched
        for coincidences.  Events falling outside this window will not be
        part of the coincidence.  Default: 200000 (i.e. 200 us).

    :return: a list of coincidences, which each consist of a list with
        indexes into the timestamps array as a pointer to the events
        making up the coincidence

    """
    coincidences = []

    # traverse all timestamps
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
            coincidences.append(c)

    return coincidences


def get_events(data, stations, coincidence, timestamps,
               get_raw_traces=False):
    """Get event data of a coincidence

    Return a list of events making up a coincidence.

    :param data: the PyTables data file
    :param stations: a list of HiSPARC event tables (normally from
        different stations, hence the name)
    :param coincidence: a coincidence, as returned by
        :func:`search_coincidences`.
    :param timestamps: the list of timestamps, as returned by
        :func:`search_coincidences`.
    :param get_raw_traces: boolean.  If true, return the compressed adc
        values instead of the uncompressed traces.

    :return: a list of tuples.  Each tuple consists of (station, event,
        traces), where event is the event row from PyTables and traces is
        a list of the uncompressed traces.

    """
    events = []

    for event in coincidence:
        timestamp, station, index = timestamps[event]
        event_table = data.getNode(stations[station], 'events')
        blob_table = data.getNode(stations[station], 'blobs')
        event = event_table[index]
        if not get_raw_traces:
            traces = get_traces(blob_table, event['traces'])
        else:
            traces = [blob_table[x] for x in event['traces']]
        events.append((stations[station], event, traces))

    return events
