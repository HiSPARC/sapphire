""" Fetch events and other data from the event summary data (ESD)

    This module enables you to access the event summary data.

    For convenience, you'll want the :func:`download_data` function.

"""
import urllib2
import urllib
import csv
import StringIO
from numpy import genfromtxt

import tables


URL = 'http://data.hisparc.nl/data/%d/events'


def download_data(file, group, station_id, start, end, func):
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
        >>> data = tables.openFile('data.h5', 'w')
        >>> sapphire.esd.download_data(data, '/s501', 501, datetime.datetime(2010, 9, 1), datetime.datetime(2010, 9, 2))

    """
    url = URL % station_id
    query_string = urllib.urlencode({'start': start, 'end': end})
    url += '?' + query_string

    data = urllib2.urlopen(url)

    func(file, group, data)


def parse1(file, group, data):
    table = create_table(file, group)
    reader = csv.reader(data, delimiter='\t')

    for line in reader:
        if line[0][0] == '#':
            continue

        (date, time, timestamp, nanoseconds, ph1, ph2, ph3, ph4, int1,
         int2, int3, int4, n1, n2, n3, n4, t1, t2, t3, t4) = line

        timestamp = int(timestamp)
        nanoseconds = int(nanoseconds)
        ext_timestamp = timestamp * int(1e9) + nanoseconds
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
        table.append([[timestamp, nanoseconds, ext_timestamp,
                       pulseheights, integrals, n1, n2, n3, n4, t1, t2,
                       t3, t4]])


def parse2(file, group, data):
    table = create_table(file, group)

    format = [('date', 'datetime64[D]'), ('time', '|S8'),
              ('timestamp', 'uint32'), ('nanoseconds', 'uint32'),
              ('pulseheights', '4int16'), ('integrals', '4int16'),
              ('n1', 'float32'), ('n2', 'float32'),
              ('n3', 'float32'), ('n4', 'float32'),
              ('t1', 'float32'), ('t2', 'float32'),
              ('t3', 'float32'), ('t4', 'float32')]


    for line in data:
        try:
            x = genfromtxt(StringIO.StringIO(line), delimiter='\t', dtype=format)
        except IOError:
            continue

        x = x.tolist()
        #table.append([[x[2], x[3], x[2] * int(1e9) + x[3], x[4], x[5],
        #               x[6], x[7], x[8], x[9], x[10], x[11], x[12],
        #               x[13]]])
        timestamp = x[2]
        nanoseconds = x[3]
        ext_timestamp = timestamp * int(1e9) + nanoseconds
        pulseheights = list(x[4])
        integrals = list(x[5])
        n1 = x[6]
        n2 = x[7]
        n3 = x[8]
        n4 = x[9]
        t1 = x[10]
        t2 = x[11]
        t3 = x[12]
        t4 = x[13]
        #print timestamp, nanoseconds, ext_timestamp, pulseheights, \
        #    integrals, n1, n2, n3, n4, t1, t2, t3, t4
        table.append([[timestamp, nanoseconds, ext_timestamp,
                       pulseheights, integrals, n1, n2, n3, n4, t1, t2,
                       t3, t4]])


def create_table(file, group):
    description = {'timestamp': tables.Time32Col(pos=0),
                   'nanoseconds': tables.UInt32Col(pos=1),
                   'ext_timestamp': tables.UInt64Col(pos=2),
                   'pulseheights': tables.Int16Col(pos=3, shape=4),
                   'integrals': tables.Int32Col(pos=4, shape=4),
                   'n1': tables.Float32Col(pos=5),
                   'n2': tables.Float32Col(pos=6),
                   'n3': tables.Float32Col(pos=7),
                   'n4': tables.Float32Col(pos=8),
                   't1': tables.Float32Col(pos=9),
                   't2': tables.Float32Col(pos=10),
                   't3': tables.Float32Col(pos=11),
                   't4': tables.Float32Col(pos=12)}

    if 'events' in file.getNode(group):
        file.removeNode(group, 'events')

    return file.createTable(group, 'events', description)
