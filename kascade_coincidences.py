""" Process HiSPARC / KASCADE coincidences

    This module reads data from the HiSPARC and KASCADE event tables and
    searches for coincidences.

"""
import datetime
import time
import os
import pylab
import operator

import tables

def do_timeshifts(h, k, shifts, dtlimit=None, limit=None):
    """Search for coincidences using multiple time shifts

    This function enables you to search for coincidences multiple times,
    using a list of time shifts. Given a data file, the events are read
    into arrays and passed on to the search_coincidences function. For
    each shift, a histogram is plotted so you can get a feel for the
    goodness of the shift. The coincidences data from the last shift is
    returned.

    :param hevents: hisparc event table
    :param kevents: kascade event table
    :param shifts: a list of time shifts
    :param dtlimit: limit on the time difference between hisparc and
        kascade events in nanoseconds.  If this limit is exceeded,
        coincidences are not stored.  Default: None.
    :param limit: an optional limit on the number of kascade events used
        in the search
    :param h: prefetched array from hisparc table (optional)
    :param k: prefetched array from kascade table (optional)

    :return: An array of coincidences from the last shift ([dt in
        nanoseconds, hisparc event id, kascade event id]).

    """
    for shift in shifts:
        print "Calculating dt's for timeshift %.9f (%d nanoseconds)" % \
              (shift, long(shift * 1e9))
        coincidences = search_coincidences_fun(h, k, shift, dtlimit)

        dt = [x[0] / 1e9 for x in coincidences]
        pylab.hist(dt, bins=100, range=(-1, 1), histtype='step',
                   label="Shift %+g s" % shift)

    finish_graph()
    return coincidences

def search_coincidences(hisparc_data, kascade_data, timeshift,
                        dtlimit=None):
    """Search for coincidences

    This function does the actual searching of coincidences. It uses a
    timeshift to shift the HiSPARC data (we know that these employ GPS
    time, so not taking UTC leap seconds into account). The shift will also
    compensate for delays in the experimental setup.

    :param hisparc_data: an array containing the hisparc data
    :param kascade_data: an array containing the kascade data
    :param timeshift: the amount of time the HiSPARC data are shifted (in
        seconds)
    :param dtlimit: limit on the time difference between hisparc and
        kascade events in nanoseconds.  If this limit is exceeded,
        coincidences are not stored.  Default: None.

    :return: An array of time differences and event ids of each KASCADE
        event and the nearest neighbour HiSPARC event.

    """
    # Shift the kascade data instead of the hisparc data. There is less of
    # it, so this is much faster.
    k = shift_data(kascade_data, -timeshift)
    h = hisparc_data

    coincidences = []

    # First loop through kascade data until we have the first event that
    # occurs _after_ the first hisparc event.
    h_idx = 0
    for k_idx in range(len(k)):
        if k[k_idx][1] > h[h_idx][1]:
            break

    while True:
        # Try to get the timestamps of the kascade event and the
        # neighbouring hisparc events.
        try:
            h_t = h[h_idx][1]
            k_t = k[k_idx][1]
            h_t_next = h[h_idx + 1][1]
        except IndexError:
            # Reached beyond the event list.
            break

        # Make sure that while the current hisparc event is _before_ the
        # kascade event, the next hisparc event should occur _after_ the
        # kascade event.  That way, the kascade event is enclosed by
        # hisparc events.
        if k_t > h_t_next:
            h_idx += 1
            continue

        # Calculate the time differences for both neighbors. Make sure to
        # get the sign right. Negative sign: the hisparc event is 'left'.
        # Positive sign: the hisparc event is 'right'.
        dt_left = h_t - k_t
        dt_right = h_t_next - k_t

        # Determine the nearest neighbor and add that to the coincidence
        # list, if dtlimit is not exceeded
        if dtlimit is None or min(abs(dt_left), abs(dt_right)) < dtlimit:
            if abs(dt_left) < abs(dt_right):
                coincidences.append((dt_left, h_idx, k_idx))
            else:
                coincidences.append((dt_right, h_idx + 1, k_idx))

        # Found a match for this kascade event, so continue with the next
        # one.
        k_idx += 1

    return coincidences

def search_coincidences_fun(hisparc_data, kascade_data, timeshift,
                            dtlimit=None):
    """Search for coincidences

    This function does the actual searching of coincidences. It uses a
    timeshift to shift the HiSPARC data (we know that these employ GPS
    time, so not taking UTC leap seconds into account). The shift will also
    compensate for delays in the experimental setup.

    :param hisparc_data: an array containing the hisparc data
    :param kascade_data: an array containing the kascade data
    :param timeshift: the amount of time the HiSPARC data are shifted (in
        seconds)
    :param dtlimit: limit on the time difference between hisparc and
        kascade events in nanoseconds.  If this limit is exceeded,
        coincidences are not stored.  Default: None.

    :return: An array of time differences and event ids of each KASCADE
        event and the nearest neighbour HiSPARC event.

    """
    # Shift the kascade data instead of the hisparc data. There is less of
    # it, so this is much faster.
    k = shift_data(kascade_data, -timeshift)
    h = hisparc_data

    coincidences = []

    # First loop through kascade data until we have the first event that
    # occurs _after_ the first hisparc event.
    h_idx = 0
    for k_idx in range(len(k)):
        if k[k_idx][1] > h[h_idx][1]:
            break

    while True:
        # Try to get the timestamps of the kascade event and the
        # neighbouring hisparc events.
        try:
            h_tt = h[h_idx - 1][1]
            h_t = h[h_idx][1]
            k_t = k[k_idx][1]
            h_t_next = h[h_idx + 1][1]
            h_tt_next = h[h_idx + 2][1]
        except IndexError:
            # Reached beyond the event list.
            break

        # Make sure that while the current hisparc event is _before_ the
        # kascade event, the next hisparc event should occur _after_ the
        # kascade event.  That way, the kascade event is enclosed by
        # hisparc events.
        if k_t > h_t_next:
            h_idx += 1
            continue

        # Calculate the time differences for both neighbors. Make sure to
        # get the sign right. Negative sign: the hisparc event is 'left'.
        # Positive sign: the hisparc event is 'right'.
        dtt_left = h_tt - k_t
        dt_left = h_t - k_t
        dt_right = h_t_next - k_t
        dtt_right = h_tt_next - k_t

        # Determine the next-to-nearest neighbor and add that to the
        # coincidence list, if dtlimit is not exceeded
        #if dtlimit is None or min(abs(dt_left), abs(dt_right)) < dtlimit:
        #    if abs(dt_left) < abs(dt_right):
        #        coincidences.append((dt_left, h_idx, k_idx))
        #    else:
        #        coincidences.append((dt_right, h_idx + 1, k_idx))
        cs = sorted([(abs(x), x) for x in (dtt_left, dt_left, dt_right,
                                   dtt_right)])
        coincidences.append((cs[1][1], cs[1][0]))

        # Found a match for this kascade event, so continue with the next
        # one.
        k_idx += 1

    return coincidences

def shift_data(data, timeshift):
    """Shift event data in time

    This function shifts the event data in time, by specifying a timeshift
    in seconds. The original data is left untouched. Returns a new array
    containing the shifted data.

    :param data: the HiSPARC or KASCADE data to be shifted
    :param timeshift: the timeshift in seconds

    :return: an array containing the original data shifted in time

    """
    # convert timeshift to an integer value in nanoseconds
    timeshift = long(timeshift * 1e9)

    return [[x[0], x[1] + timeshift] for x in data]

def finish_graph():
    """Finish the histogram

    This function places a legend, axes titles and the like on the current
    figure.

    """
    pylab.legend()
    pylab.xlabel("Time difference (s)")
    pylab.ylabel("Counts")
    pylab.title("Nearest neighbour events for HiSPARC / KASCADE")
    pylab.gca().axis('auto')
    pylab.gcf().show()


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'r')

    try:
        h, k
    except NameError:
        h = data.root.datasets.h.read()
        k = data.root.datasets.knew.read()

    #do_timeshifts(h, k, [-12, -13, -14], limit=1000)
    do_timeshifts(h, k, [0, 1, 2])
    c = search_coincidences_fun(h, k, 0)
