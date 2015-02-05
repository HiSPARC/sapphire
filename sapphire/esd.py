""" Fetch events and other data from the event summary data (ESD)

    This module enables you to access the event summary data.

    If you are in a real hurry and know what you're doing (and took the
    time to read this far), you can call the :func:`quick_download`
    function like this::

        >>> import sapphire.esd
        >>> data = sapphire.esd.quick_download(501)

    For regular use, look up :func:`download_data`.

"""
import urllib2
import urllib
import csv
import os.path
import calendar
import time
import datetime
import itertools
import collections

import tables
from progressbar import ProgressBar, ETA, Bar, Percentage

from . import api
from . import clusters
from . import storage


EVENTS_URL = 'http://data.hisparc.nl/data/{station_number:d}/events?{query}'
WEATHER_URL = 'http://data.hisparc.nl/data/{station_number:d}/weather?{query}'
COINCIDENCES_URL = 'http://data.hisparc.nl/data/network/coincidences/?{query}'


def quick_download(station_number, date=None):
    """Quickly download some data

    :param station_number: The HiSPARC station number
    :param date: the date for which to get data (datetime.datetime
        instance), passed unchanged to the :func:`download_data` as
        :param start:.
    :returns: handle to an open PyTables file

    Everything is handled by this function, including file creation.
    Expect no frills: you just get yesterday's data.

    Example usage::

        >>> import sapphire.esd
        >>> data = sapphire.esd.quick_download(501)
        >>> print data
        data1.h5 (File) u''
        Last modif.: 'Mon Jun  9 22:03:50 2014'
        Object Tree:
        / (RootGroup) u''
        /s501 (Group) u''
        /s501/events (Table(58898,)) ''

    """
    path = _first_available_numbered_path()
    data = tables.open_file(path, 'w')
    download_data(data, None, station_number, date)
    return data


def _first_available_numbered_path():
    """Find first available file name in sequence

    If data1.h5 is taken, return data2.h5, etc.

    """
    path = 'data%d.h5'
    return next(path % idx for idx in itertools.count(start=1)
                if not os.path.exists(path % idx))


def load_data(file, group, csv_file, type='events'):
    """Download event summary data

    :param file: the PyTables datafile handler
    :param group: the PyTables destination group, which need not exist
    :param csv_file: path to the csv file downloaded from the HiSPARC
                     Public Database
    :param type: the datatype to download, either 'events' or 'weather'

    Example::

        >>> import tables
        >>> import sapphire.esd
        >>> data = tables.open_file('data.h5', 'w')
        >>> sapphire.esd.load_data(data, '/s501', 'events-s501-20130910.csv')

    """
    if type == 'events':
        table = _get_or_create_events_table(file, group)
        read_and_store = _read_line_and_store_event
    elif type == 'weather':
        table = _get_or_create_weather_table(file, group)
        read_and_store = _read_line_and_store_weather
    else:
        raise ValueError("Data type not recognized.")

    with open(csv_file, 'rb') as data:
        reader = csv.reader(data, delimiter='\t')
        for line in reader:
            read_and_store(line, table)
        table.flush()


def download_data(file, group, station_number, start=None, end=None,
                  type='events'):
    """Download event summary data

    :param file: the PyTables datafile handler
    :param group: the PyTables destination group, which need not exist
    :param station_number: The HiSPARC station number for which to get data
    :param start: a datetime instance defining the start of the search
        interval
    :param end: a datetime instance defining the end of the search
        interval
    :param type: the datatype to download, either 'events' or 'weather'.

    If group is None, use '/s<station_number>' as a default.

    The start and stop parameters may both be None.  In that case,
    yesterday's data is downloaded.  If only end is None, a single day's
    worth of data is downloaded, starting at the datetime specified with
    start.

    Example::

        >>> import tables
        >>> import datetime
        >>> import sapphire.esd
        >>> data = tables.open_file('data.h5', 'w')
        >>> sapphire.esd.download_data(data, '/s501', 501,
        ...     datetime.datetime(2013, 9, 1), datetime.datetime(2013, 9, 2))

    """
    # sensible default for group name
    if group is None:
        group = '/s%d' % station_number

    # sensible defaults for start and end
    if start is None:
        if end is not None:
            raise RuntimeError("Start is None, but end is not. "
                               "I can't go on like this.")
        else:
            yesterday = datetime.date.today() - datetime.timedelta(days=1)
            start = datetime.datetime.combine(yesterday, datetime.time(0, 0))
    if end is None:
        end = start + datetime.timedelta(days=1)

    # build and open url, create tables and set read function
    query = urllib.urlencode({'start': start, 'end': end})
    if type == 'events':
        url = EVENTS_URL.format(station_number=station_number, query=query)
        table = _get_or_create_events_table(file, group)
        read_and_store = _read_line_and_store_event
    elif type == 'weather':
        url = WEATHER_URL.format(station_number=station_number, query=query)
        table = _get_or_create_weather_table(file, group)
        read_and_store = _read_line_and_store_weather
    else:
        raise ValueError("Data type not recognized.")

    data = urllib2.urlopen(url)

    # keep track of event timestamp within [start, end] interval for
    # progressbar
    t_start = calendar.timegm(start.utctimetuple())
    t_end = calendar.timegm(end.utctimetuple())
    t_delta = t_end - t_start
    pbar = ProgressBar(maxval=1., widgets=[Percentage(), Bar(), ETA()]).start()

    # loop over lines in csv as they come streaming in
    prev_update = time.time()
    reader = csv.reader(data, delimiter='\t')
    for line in reader:
        timestamp = read_and_store(line, table)
        # update progressbar every .5 seconds
        if time.time() - prev_update > .5 and not timestamp == 0.:
            pbar.update((1. * timestamp - t_start) / t_delta)
            prev_update = time.time()
    table.flush()
    pbar.finish()


def download_coincidences(file, cluster=None, stations=None,
                          start=None, end=None, n=2):
    """Download event summary data coincidences

    :param file: The PyTables datafile handler.
    :param cluster: The HiSPARC cluster name for which to get data.
    :param stations: A list of HiSPARC station numbers for which to get data.
    :param start: a datetime instance defining the start of the search
        interval.
    :param end: a datetime instance defining the end of the search
        interval.
    :param n: the minimum number of events in the coincidence.

    The start and end parameters may both be None.  In that case,
    yesterday's data is downloaded.  If only end is None, a single day's
    worth of data is downloaded, starting at the datetime specified with
    start.

    Optionally either a cluster or stations can be defined to limit the
    results to include only events from those stations.

    Example::

        import tables
        import datetime
        import sapphire.esd
        data = tables.open_file('data_coincidences.h5', 'w')
        sapphire.esd.download_coincidences(data, cluster='Aarhus',
            start=datetime.datetime(2013, 9, 1),
            end=datetime.datetime(2013, 9, 2), n=3)

    """
    # sensible defaults for start and end
    if start is None:
        if end is not None:
            raise RuntimeError("Start is None, but end is not. "
                               "I can't go on like this.")
        else:
            yesterday = datetime.date.today() - datetime.timedelta(days=1)
            start = datetime.datetime.combine(yesterday, datetime.time(0, 0))
    if end is None:
        end = start + datetime.timedelta(days=1)

    # build and open url, create tables and set read function
    query = urllib.urlencode({'cluster': cluster, 'stations': stations,
                              'start': start, 'end': end, 'n': n})
    url = COINCIDENCES_URL.format(query=query)
    station_groups = _get_station_groups()
    table = _get_or_create_coincidences_tables(file, station_groups)
    station_numbers = _get_or_create_station_numbers(table)

    data = urllib2.urlopen(url, timeout=1800)

    # keep track of event timestamp within [start, end] interval for
    # progressbar
    t_start = calendar.timegm(start.utctimetuple())
    t_end = calendar.timegm(end.utctimetuple())
    t_delta = t_end - t_start
    pbar = ProgressBar(maxval=1., widgets=[Percentage(), Bar(), ETA()]).start()

    # loop over lines in csv as they come streaming in, keep temporary
    # lists untill a full coincidence is in.
    prev_update = time.time()
    reader = csv.reader(data, delimiter='\t')
    current_coincidence = 0
    coincidence = []
    for line in reader:
        if line[0][0] == '#':
            continue
        elif int(line[0]) == current_coincidence:
            station_numbers.add(int(line[1]))
            coincidence.append(line)
        else:
            station_numbers.add(int(line[1]))
            # Full coincidence has been received, store it.
            timestamp = _read_lines_and_store_coincidence(file,
                                                          coincidence,
                                                          station_groups)
            # update progressbar every .5 seconds
            if time.time() - prev_update > .5 and not timestamp == 0.:
                pbar.update((1. * timestamp - t_start) / t_delta)
                prev_update = time.time()
            coincidence = [line]
            current_coincidence = int(line[0])
            file.flush()

    if len(coincidence):
        # Store last coincidence
        _read_lines_and_store_coincidence(file, coincidence, station_groups)

    pbar.finish()

    cluster = clusters.HiSPARCStations(station_numbers)
    table._v_attrs.cluster = cluster
    file.flush()


def _get_or_create_station_numbers(table):
    """Get station numbers from existing cluster attribute or a new set

    :param table: coincidence table in PyTables file.
    :returns: set including existing stations in a cluster.

    """
    try:
        cluster = table._v_attrs.cluster
        return {s.number for s in cluster.stations}
    except AttributeError:
        return set()


def _get_station_groups():
    """Generate groups names for all stations

    Use the same hierarchy (cluster/station) as used in the HiSPARC
    datastore.

    """
    groups = collections.OrderedDict()
    network = api.Network()
    clusters = network.clusters()
    s_index = 0
    for cluster in clusters:
        stations = api.Network().station_numbers(cluster=cluster['number'])
        for station in stations:
            groups[station] = {'group': ('/hisparc/cluster_%s/station_%d' %
                                         (cluster['name'].lower(), station)),
                               's_index': s_index}
            s_index += 1
    return groups


def _get_or_create_coincidences_tables(file, station_groups):
    """Get or create event table in PyTables file"""

    try:
        return file.get_node('/', 'coincidences')
    except tables.NoSuchNodeError:
        return _create_coincidences_tables(file, station_groups)


def _create_coincidences_tables(file, station_groups):
    """Setup coincidence tables"""

    group = '/coincidences'

    # Create coincidences table
    description = storage.Coincidence
    s_columns = {'s%d' % station: tables.BoolCol(pos=p)
                 for p, station in enumerate(station_groups.iterkeys(), 12)}
    description.columns.update(s_columns)
    coincidences = file.create_table(group, 'coincidences', description,
                                     createparents=True)

    # Create c_index
    file.create_vlarray(group, 'c_index', tables.UInt32Col(shape=2))

    # Create and fill s_index
    s_index = file.create_vlarray(group, 's_index', tables.VLStringAtom())
    for station_group in station_groups.itervalues():
        s_index.append(station_group['group'])

    return coincidences


def _get_or_create_events_table(file, group):
    """Get or create event table in PyTables file"""

    try:
        return file.get_node(group, 'events')
    except tables.NoSuchNodeError:
        return _create_events_table(file, group)


def _create_events_table(file, group):
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

    return file.create_table(group, 'events', description, createparents=True)


def _get_or_create_weather_table(file, group):
    """Get or create event table in PyTables file"""

    try:
        return file.get_node(group, 'weather')
    except tables.NoSuchNodeError:
        return _create_weather_table(file, group)


def _create_weather_table(file, group):
    """Create weather table in PyTables file

    Create a weather table containing the ESD weather columns which are
    available in the CSV download.

    :param file: PyTables file
    :param group: the group to contain the weather table, which need not
                  exist

    """
    description = {'event_id': tables.UInt32Col(pos=0),
                   'timestamp': tables.Time32Col(pos=1),
                   'temp_inside': tables.Float32Col(pos=2),
                   'temp_outside': tables.Float32Col(pos=3),
                   'humidity_inside': tables.Int16Col(pos=4),
                   'humidity_outside': tables.Int16Col(pos=5),
                   'barometer': tables.Float32Col(pos=6),
                   'wind_dir': tables.Int16Col(pos=7),
                   'wind_speed': tables.Int16Col(pos=8),
                   'solar_rad': tables.Int16Col(pos=9),
                   'uv': tables.Int16Col(pos=10),
                   'evapotranspiration': tables.Float32Col(pos=11),
                   'rain_rate': tables.Float32Col(pos=12),
                   'heat_index': tables.Int16Col(pos=13),
                   'dew_point': tables.Float32Col(pos=14),
                   'wind_chill': tables.Float32Col(pos=15)}

    return file.create_table(group, 'weather', description, createparents=True)


def _read_lines_and_store_coincidence(file, coincidence, station_groups):
    """Read CSV lines and store coincidence

    Read lines from the CSV download and store the coincidence and events.
    Return the coincidence timestamp to keep track of the progress.

    :param coincidence: text lines from the CSV file for one coincidence
    :param file: pytables file for storage
    :return: coincidence timestamp

    """
    c_idx = []
    coincidences = file.get_node('/coincidences', 'coincidences')
    row = coincidences.row
    row['id'] = len(coincidences)
    row['N'] = len(coincidence)
    row['timestamp'] = int(coincidence[0][4])
    row['nanoseconds'] = int(coincidence[0][5])
    row['ext_timestamp'] = (int(coincidence[0][4]) * int(1e9) +
                            int(coincidence[0][5]))
    for event in coincidence:
        station_number = int(event[1])
        row['s%d' % station_number] = True
        group_path = station_groups[station_number]['group']
        group = _get_or_create_events_table(file, group_path)
        s_idx = station_groups[station_number]['s_index']
        e_idx = len(group)
        c_idx.append((s_idx, e_idx))
        _read_line_and_store_event(event[2:], group)

    row.append()
    c_index = file.get_node('/coincidences', 'c_index')

    c_index.append(c_idx)

    return int(coincidence[0][4])


def _read_line_and_store_event(line, table):
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

    row = table.row

    # convert string values to correct data types or calculate values
    row['event_id'] = len(table)
    row['timestamp'] = int(timestamp)
    row['nanoseconds'] = int(nanoseconds)
    row['ext_timestamp'] = int(timestamp) * int(1e9) + int(nanoseconds)
    row['pulseheights'] = [int(ph1), int(ph2), int(ph3), int(ph4)]
    row['integrals'] = [int(int1), int(int2), int(int3), int(int4)]
    row['n1'] = float(n1)
    row['n2'] = float(n2)
    row['n3'] = float(n3)
    row['n4'] = float(n4)
    row['t1'] = float(t1)
    row['t2'] = float(t2)
    row['t3'] = float(t3)
    row['t4'] = float(t4)
    row['t_trigger'] = float(t_trigger)

    # store event
    row.append()

    return int(timestamp)


def _read_line_and_store_weather(line, table):
    """Read CSV line and store weather data

    Read a line from the CSV download and store weather.  Return the
    weather timestamp to keep track of the progress.

    :param line: text line from the CSV file
    :param table: pytables table for weather storage
    :return: weather timestamp

    """
    # ignore comment lines
    if line[0][0] == '#':
        return 0.

    # break up CSV line
    (date, time, timestamp, temperature_inside, temperature_outside,
     humidity_inside, humidity_outside, atmospheric_pressure,
     wind_direction, wind_speed, solar_radiation, uv_index,
     evapotranspiration, rain_rate, heat_index, dew_point, wind_chill) = line

    row = table.row

    # convert string values to correct data types
    row['event_id'] = len(table)
    row['timestamp'] = int(timestamp)
    row['temp_inside'] = float(temperature_inside)
    row['temp_outside'] = float(temperature_outside)
    row['humidity_inside'] = int(humidity_inside)
    row['humidity_outside'] = int(humidity_outside)
    row['barometer'] = float(atmospheric_pressure)
    row['wind_dir'] = int(wind_direction)
    row['wind_speed'] = int(wind_speed)
    row['solar_rad'] = int(solar_radiation)
    row['uv'] = int(uv_index)
    row['evapotranspiration'] = float(evapotranspiration)
    row['rain_rate'] = float(rain_rate)
    row['heat_index'] = int(heat_index)
    row['dew_point'] = float(dew_point)
    row['wind_chill'] = float(wind_chill)

    # store event
    row.append()

    return int(timestamp)
