""" Search for coincidences between HiSPARC stations

    This module can be used to search for coincidences between several
    HiSPARC stations.

    To search for coincidences, use the :func:`search_coincidences`
    function.
"""

import numpy as np
import time

COINC_WINDOW = long(200e-6 * 1e9)

def search_coincidences(data, stations, shifts=None, limit=None):
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
    :param shifts: a list of time shifts which may contain 'None'.
    :param limit: limit the number of events which are processed

    :return: coincidences, timestamps. First a list of coincidences, which
        each consist of a list with indexes into the timestamps array as a
        pointer to the events making up the coincidence. Then, a list of
        tuples.  Each tuple consists of a timestamp followed by an index
        into the stations list which designates the detector
        station which measured the event.

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
        ((1235433610410730837, 0), (1235433610410731004, 2))

    """
    # get the 'events' tables from the groups or groupnames
    stations = [data.getNode(x, 'events') for x in stations]

    # calculate the shifts in nanoseconds and cast them to long.
    # (prevent upcasting timestamps to float64 further on)
    if shifts:
        for i in range(len(shifts)):
            if shifts[i]:
                shifts[i] = long(shifts[i] * 1e9)

    timestamps = retrieve_timestamps(stations, shifts, limit)
    coincidences = do_search_coincidences(timestamps)

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
        station which measured the event.

    """
    timestamps = []
    for i in range(len(stations)):
        ts = [(x['ext_timestamp'], i) for x in stations[i][:limit]]
        try:
            # shift data. carefully avoid upcasting (we're adding two
            # longs, which is a long, and casting that back to uint64. if
            # we're not careful, an intermediate value will be a float64,
            # which doesn't hold the precision to store nanoseconds.
            ts = [(np.uint64(long(x[0]) + shifts[i]), x[1]) for x in ts]
        except (TypeError, IndexError):
            # shift is None or doesn't exist
            pass
        timestamps.extend(ts)

    # sort the timestamps
    timestamps.sort()

    return timestamps

def do_search_coincidences(timestamps):
    """Search for coincidences in a set of timestamps

    Given a set of timestamps, search for coincidences.  That is, search
    for events which occured almost at the same time and thus might be the
    result of an extended air shower.

    :param timestamps: a list of tuples (timestamp, station_idx) which
        will be searched

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
            if timestamps[j][0] - t0 < COINC_WINDOW:
                c.append(j)
            else:
                # coincidence window has passed, break for-loop
                break

        # if we have more than one event in the coincidence, save it
        if len(c) > 1:
            coincidences.append(c)

    return coincidences
