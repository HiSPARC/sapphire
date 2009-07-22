""" Process HiSPARC / KASCADE coincidences

    This module reads data from the HiSPARC and KASCADE event tables and
    searches for coincidences.

"""
import datetime
import time
import os

def do_timeshifts(datafile, shifts, limit=None):
    """Search for coincidences using multiple time shifts

    This function enables you to search for coincidences multiple times,
    using a list of time shifts. Given a data file, the events are read
    into arrays and passed on to the search_coincidences function. For
    each shift, a histogram is plotted so you can get a feel for the
    goodness of the shift. The coincidences data from the last shift is
    returned.

    Arguments:
    datafile    the data file containing the events
    shifts      a list of time shifts
    limit       an optional limit on the number of kascade events used in
                the search

    Returns:
    An array of coincidences from the last shift ([dt in nanoseconds,
    hisparc event id, kascade event id]).

    """
    # Get arrays from the tables. This is much, much faster than working
    # from the tables directly. Pity.
    h, k = get_arrays_from_tables(datafile.root.hisparc.events,
                                  datafile.root.kascade.events, limit)

    for shift in shifts:
        print "Calculating dt's for timeshift %.9f (%d nanoseconds)" % \
              (shift, long(shift * 1e9))
        coincidences = search_coincidences(h, k, shift)

        dt = [x[0] / 1e9 for x in coincidences]
        hist(dt, bins=100, range=(-1, 1), histtype='step',
             label="Shift %+g s" % shift)

    finish_graph()
    return coincidences

def store_coincidences(datafile, coincidences):
    """Store coincidences in a table

    This function stores coincidences which are found by
    search_coincidences in a table, so data can be easily retrieved without
    resorting to lookups which span multiple tables.

    Arguments:
    datafile            datafile to hold the coincidences
    coincidences        a list of coincidences, as given by
                        search_coincidences

    """
    table = datafile.root.coincidences.events
    old_data_length = len(table)
    tablerow = table.row

    for coincidence in coincidences:
        hisparc = datafile.root.hisparc.events[coincidence[1]]
        kascade = datafile.root.kascade.events[coincidence[2]]
        tablerow['hisparc_event_id'] = hisparc['event_id']
        tablerow['kascade_event_id'] = kascade['event_id']
        tablerow['hisparc_timestamp'] = hisparc['timestamp']
        tablerow['hisparc_nanoseconds'] = hisparc['nanoseconds']
        tablerow['hisparc_ext_timestamp'] = hisparc['ext_timestamp']
        tablerow['hisparc_pulseheights'] = hisparc['pulseheights']
        tablerow['hisparc_integrals'] = hisparc['integrals']
        tablerow['kascade_timestamp'] = kascade['timestamp']
        tablerow['kascade_nanoseconds'] = kascade['nanoseconds']
        tablerow['kascade_ext_timestamp'] = kascade['ext_timestamp']
        tablerow['kascade_energy'] = kascade['energy']
        tablerow['kascade_core_pos'] = kascade['core_pos']
        tablerow['kascade_zenith'] = kascade['zenith']
        tablerow['kascade_azimuth'] = kascade['azimuth']
        tablerow['kascade_Num_e'] = kascade['Num_e']
        tablerow['kascade_Num_mu'] = kascade['Num_mu']
        tablerow['kascade_dens_e'] = kascade['dens_e']
        tablerow['kascade_dens_mu'] = kascade['dens_mu']
        tablerow['kascade_P200'] = kascade['P200']
        tablerow['kascade_T200'] = kascade['T200']
        tablerow.append()
    table.flush()

    # Flush old data
    if old_data_length:
        print "Flushing old data..."
        table.removeRows(0, old_data_length)

def search_coincidences(hisparc_data, kascade_data, timeshift, limit=None):
    """Search for coincidences

    This function does the actual searching of coincidences. It uses a
    timeshift to shift the HiSPARC data (we know that these employ GPS
    time, so not taking UTC leap seconds into account). The shift will also
    compensate for delays in the experimental setup.

    Arguments:
    hisparc_data        an array containing the hisparc data
    kascade_data        an array containing the kascade data
    timeshift           the amount of time the HiSPARC data are shifted (in
                        seconds)

    Output:
    An array of time differences and event ids of each KASCADE event and
    the nearest neighbour HiSPARC event.

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
        # list.
        if abs(dt_left) < abs(dt_right):
            coincidences.append((dt_left, h_idx, k_idx))
        else:
            coincidences.append((dt_right, h_idx + 1, k_idx))

        # Found a match for this kascade event, so continue with the next
        # one.
        k_idx += 1

    return coincidences

def shift_data(data, timeshift):
    """Shift event data in time

    This function shifts the event data in time, by specifying a timeshift
    in seconds. The original data is left untouched. Returns a new array
    containing the shifted data.

    Arguments:
    data        the HiSPARC or KASCADE data to be shifted
    timeshift   the timeshift in seconds

    Returns:
    an array containing the original data shifted in time

    """
    # convert timeshift to an integer value in nanoseconds
    timeshift = long(timeshift * 1e9)

    return [[x[0], x[1] + timeshift] for x in data]

def finish_graph():
    """Finish the histogram

    This function places a legend, axes titles and the like on the current
    figure.

    """
    legend()
    xlabel("Time difference (s)")
    ylabel("Counts")
    title("Nearest neighbour events for HiSPARC / KASCADE")
    gca().axis('auto')
    gcf().show()

def get_arrays_from_tables(h, k, limit):
    """Get data arrays from data tables

    This function returns an array of values extracted from the event
    tables with hisparc and kascade data. It honors a limit and only
    fetches events which fall inside the time window.

    Caveat: because the timeshift is not yet known, a few coincidences
    may fall outside the time window and not be taken into account.

    Arguments:
    h           hisparc event table
    k           kascade event table
    limit       limit on the number of kascade events

    Returns:
    Two arrays containing hisparc and kascade data ([event id, timestamp in
    nanoseconds])

    """
    try:
        k_t = k[limit - 1]['timestamp']
    except (IndexError, TypeError):
        k_t = k[-1]['timestamp']
    h_t = h[-1]['timestamp']

    t_end = min([k_t, h_t])

    k = [[x['event_id'], x['ext_timestamp']] for x in \
         k.where('timestamp <= t_end')]
    h = [[x['event_id'], x['ext_timestamp']] for x in \
         h.where('timestamp <= t_end')]

    return h, k


if __name__ == '__main__':
    print "Careful: the following search is limited to 1000 kascade events"
    print "The complete statement would be:"
    print "c = do_timeshifts(data, [-13.180213654])"
    c = do_timeshifts(data, [-13.180213654], limit=1000)
