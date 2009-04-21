""" Process HiSPARC / KASCADE coincidences

    This module reads data from the HiSPARC and KASCADE event tables and
    searches for coincidences.

"""
import datetime
import time
import os

def search_coincidences(datafile, timeshift, limit=None, store=False):
    """Search for coincidences and optionally store the results

    This function does the actual searching of coincidences. It uses a
    timeshift to shift the HiSPARC data (we know that these employ GPS
    time, so not taking UTC leap seconds into account). The shift will also
    compensate for delays in the experimental setup.

    Storing the results is only optional, to be able to play with different
    timeshifts. To make this possible, an array of all the time differences
    is returned to assess the goodness of the timeshift.

    Arguments:
    datafile    an instance of a pytables data file
    timeshift   the amount of time the HiSPARC data are shifted (in
                seconds)
    limit       if given, the maximum number of kascade events used in the
                search
    store       a boolean to control whether or not the results are stored
                in the data file

    Output:
    An array of time differences between each KASCADE event and the nearest
    neighbour HiSPARC event.

    """
    h = data.root.hisparc.events
    k = data.root.kascade.events
    # If limit is given, use only that number of kascade events
    if limit:
        k = k[:limit]

    # Convert the timeshift to nanoseconds
    timeshift *= 1e9

    coincidences = []
    dt_list = []

    # First loop through kascade data until we have the first event that
    # occurs _after_ the first hisparc event.
    h_idx = 0
    for k_idx in range(len(k)):
        if k[k_idx]['timestamp'] > h[h_idx]['timestamp']:
            break

    while True:
        # Try to get the timestamps of the kascade event and the
        # neighbouring hisparc events. Calculate the timestamps in
        # nanoseconds, for maximum precision
        try:
            h_t = h[h_idx]['timestamp'] * 1e9 + h[h_idx]['nanoseconds'] + \
                  timeshift
            k_t = k[k_idx]['timestamp'] * 1e9 + k[k_idx]['nanoseconds']
            h_t_next = h[h_idx + 1]['timestamp'] * 1e9 + \
                       h[h_idx + 1]['nanoseconds'] + timeshift
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
        # list and the time differences list.
        if abs(dt_left) < abs(dt_right):
            dt_list.append(dt_left)
            coincidences.append((dt_left, h_idx, k_idx))
        else:
            dt_list.append(dt_right)
            coincidences.append((dt_right, h_idx + 1, k_idx))

        # Found a match for this kascade event, so continue with the next
        # one.
        k_idx += 1

    if store:
        print "Not implemented yet!"

    return dt_list

def do_timeshifts(datafile, shifts, limit=None):
    for shift in shifts:
        print "Calculating dt's for timeshift", shift
        dt = search_coincidences(datafile, shift, limit, store=False)
        dt = [x / 1e9 for x in dt]
        hist(dt, bins=100, range=(-1, 1), histtype='step',
             label="Shift %+g s" % shift)
    finish_graph()
    return dt

def finish_graph():
    legend()
    xlabel("Time difference (s)")
    ylabel("Counts")
    title("Nearest neighbour events for HiSPARC / KASCADE")
    gca().axis('auto')
    gcf().show()

#if __name__ == '__main__':
#    do_timeshifts(data, [13.180212926864623])
