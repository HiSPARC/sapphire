""" Fetch events and other data from the event summary data (ESD)

    This module enables you to access the event summary data.

    For convenience, you'll want the :func:`download_data` function.

"""
import urllib2
import urllib
import csv
import os.path
import calendar
import time

import tables
import progressbar


URL = 'http://data.hisparc.nl/data/%d/events'


def download_data(file, group, station_id, start, end):
    """Download event summary data

    :param file: The PyTables datafile handler
    :param group: The PyTables destination group, which need not exist
    :param station_id: The HiSPARC station number for which to get events
    :param start: a datetime instance defining the start of the search
        interval
    :param end: a datetime instance defining the end of the search
        interval

    Example::

        >>> import tables
        >>> import datetime
        >>> import sapphire.esd
        >>> data = tables.open_file('data.h5', 'w')
        >>> sapphire.esd.download_data(data, '/s501', 501,
        ... datetime.datetime(2013, 9, 1), datetime.datetime(2013, 9, 2))

    """
    # build and open url
    url = URL % station_id
    query_string = urllib.urlencode({'start': start, 'end': end})
    url += '?' + query_string
    data = urllib2.urlopen(url)

    # keep track of event timestamp within [start, end] interval for
    # progressbar
    t_start = calendar.timegm(start.utctimetuple())
    t_end = calendar.timegm(end.utctimetuple())
    t_delta = t_end - t_start
    pbar = progressbar.ProgressBar(maxval=1.,
                                   widgets=[progressbar.Percentage(),
                                            progressbar.Bar(),
                                            progressbar.ETA()]).start()

    # create events table
    table = create_table(file, group)

    # event loop
    prev_update = time.time()
    reader = csv.reader(data, delimiter='\t')
    for line in reader:
        timestamp = read_line_and_store_event(line, table)

        # update progressbar every .5 seconds
        if time.time() - prev_update > .5 and not timestamp == 0.:
            pbar.update((1. * timestamp - t_start) / t_delta)
            prev_update = time.time()
    pbar.finish()


def create_table(file, group):
    """Create event table in PyTables file

    Create an event table containing the ESD data columns which are
    available in the CSV download.

    :param file: PyTables file
    :param group: the group to contain the events table, which need not
        exist

    """
    description = {'event_id': tables.UInt32Col(pos=0),
                   'timestamp': tables.Time32Col(pos=1),
                   'nanoseconds': tables.UInt32Col(pos=2),
                   'ext_timestamp': tables.UInt64Col(pos=3),
                   'pulseheights': tables.Int16Col(pos=4, shape=4),
                   'integrals': tables.Int32Col(pos=5, shape=4),
                   'n1': tables.Float32Col(pos=6),
                   'n2': tables.Float32Col(pos=7),
                   'n3': tables.Float32Col(pos=8),
                   'n4': tables.Float32Col(pos=9),
                   't1': tables.Float32Col(pos=10),
                   't2': tables.Float32Col(pos=11),
                   't3': tables.Float32Col(pos=12),
                   't4': tables.Float32Col(pos=13),
                   't_trigger': tables.Float32Col(pos=14)}

    if group not in file:
        head, tail = os.path.split(group)
        file.create_group(head, tail)

    return file.create_table(group, 'events', description)


def read_line_and_store_event(line, table):
    """Read CSV line and store event

    Read a line from the CSV download and store event.  Return the event
    timestamp to keep track of the progress.

    :param line: text line from the CSV file
    :param table: pytables table for event storage
    :return: event timestamp

    """
    # ignore comment lines
    if line[0][0] == '#':
        return 0.

    # break up CSV line
    (date, time_str, timestamp, nanoseconds, ph1, ph2, ph3, ph4, int1,
     int2, int3, int4, n1, n2, n3, n4, t1, t2, t3, t4, t_trigger) = line

    # convert string values to correct data types or calculate values
    event_id = len(table)
    timestamp = int(timestamp)
    nanoseconds = int(nanoseconds)
    ext_timestamp = timestamp * 1000000000 + nanoseconds
    pulseheights = [int(ph1), int(ph2), int(ph3), int(ph4)]
    integrals = [int(int1), int(int2), int(int3), int(int4)]
    n1 = float(n1)
    n2 = float(n2)
    n3 = float(n3)
    n4 = float(n4)
    t1 = float(t1)
    t2 = float(t2)
    t3 = float(t3)
    t4 = float(t4)
    t_trigger = float(t_trigger)

    # store event
    table.append([[event_id, timestamp, nanoseconds, ext_timestamp,
                   pulseheights, integrals, n1, n2, n3, n4, t1, t2,
                   t3, t4, t_trigger]])

    return timestamp
