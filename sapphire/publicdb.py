""" Fetch events and other data from the public database

    This module enables you to access the public database and even the raw
    event data.

    For convenience, you'll want the :func:`download_data` function.

"""
import xmlrpclib
import urllib
import datetime
import tables
import os
import calendar
import re

import logging

logger = logging.getLogger('hisparc.publicdb')

#PUBLICDB_XMLRPC_URL = 'http://localhost:8000/raw_data/rpc'
PUBLICDB_XMLRPC_URL = 'http://data.hisparc.nl/raw_data/rpc'


def download_data(file, group, station_id, start, end, get_blobs=False):
    """Download raw data from the datastore

    This function downloads data from the datastore, using the XML-RPC API
    exposed by the public database.

    :param file: The PyTables datafile handler
    :param group: The PyTables destination group, which need not exist
    :param station_id: The HiSPARC station number for which to get events
    :param start: a datetime instance defining the start of the search
        interval
    :param end: a datetime instance defining the end of the search
        interval
    :param get_blobs: boolean, select whether binary data like traces
        should be fetched

    Example::

        >>> import tables
        >>> import datetime
        >>> import sapphire.publicdb
        >>> data = tables.open_file('data.h5', 'w')
        >>> start = datetime.datetime(2010, 9, 1)
        >>> end = datetime.datetime(2010, 9, 2)
        >>> sapphire.publicdb.download_data(data, '/s501', 501, start, end)
        INFO:hisparc.publicdb:2010-09-01 00:00:00 None
        INFO:hisparc.publicdb:Getting server data URL (2010-09-01 00:00:00)
        INFO:hisparc.publicdb:Downloading data...
        INFO:hisparc.publicdb:Storing data...
        INFO:hisparc.publicdb:Done.

    """
    server = xmlrpclib.ServerProxy(PUBLICDB_XMLRPC_URL)

    for t0, t1 in datetimerange(start, end):
        logger.info("%s %s" % (t0, t1))
        logger.info("Getting server data URL (%s)" % t0)
        try:
            url = server.hisparc.get_data_url(station_id, t0, get_blobs)
        except Exception, exc:
            if re.search("No data", str(exc)):
                logger.warning("No data for %s" % t0)
                continue
            else:
                raise
        logger.info("Downloading data...")
        tmp_datafile, headers = urllib.urlretrieve(url)
        logger.info("Storing data...")
        _store_data(file, group, tmp_datafile, t0, t1)
        logger.info("Done.")


def _store_data(dst_file, dst_group, src_filename, t0, t1):
    """Copy data from a temporary file to the destination file

    This function takes a file containing downloaded data and copies it to
    the destination file, based on start and end timestamps.

    This can be further optimized at the expense of spaghetti code.

    """
    # Open in rw mode, need to update blob idxs, if necessary
    with tables.open_file(src_filename, 'a') as src_file:
        src_group = src_file.list_nodes('/')[0]
        dst_group = _get_or_create_group(dst_file, dst_group)

        if 'blobs' in dst_group:
            len_blobs = len(dst_file.get_node(dst_group, 'blobs'))
        else:
            len_blobs = 0

        for node in src_file.list_nodes(src_group):
            dst_node = _get_or_create_node(dst_file, dst_group, node)

            if node.name == 'blobs':
                for row in node:
                    dst_node.append(row)

            elif node.name == 'events':
                if not t1:
                    cond = 'timestamp >= %d' % \
                        calendar.timegm(t0.utctimetuple())
                else:
                    cond = '(%d <= timestamp) & (timestamp <= %d)' % \
                        (calendar.timegm(t0.utctimetuple()),
                         calendar.timegm(t1.utctimetuple()))

                if len_blobs:
                    for row in node.read_where(cond):
                        row['traces'] += len_blobs
                        dst_node.append([tuple(row)])
                else:
                    for row in node.read_where(cond):
                        dst_node.append([tuple(row)])

            elif node.name == 'errors' and len_blobs:
                for row in node:
                    row['messages'] += len_blobs
                    dst_node.append([row[:]])

            elif node.name == 'config' and len_blobs:
                for row in node:
                    row['mas_version'] += len_blobs
                    row['slv_version'] += len_blobs
                    row['password'] += len_blobs
                    row['buffer'] += len_blobs
                    dst_node.append([row[:]])

            else:
                rows = node.read()
                dst_node.append(rows)

    os.remove(src_filename)
    dst_file.flush()


def datetimerange(start, stop):
    """Generator for datetime ranges

    This is a very specific generator for datetime ranges.  Based on a
    start and stop value, it generates one day intervals.  On the first
    and subsequent full days, it only returns a start value.  On the last
    day, it returns a start and stop value.  See example below.

    :param start: a datetime instance
    :param stop: a datetime instance

    Example::

        >>> for x in datetimerange(datetime.datetime(2010, 1, 1, 11),
        ...                        datetime.datetime(2010, 1, 1, 13)):
        ...     x
        ...
        (datetime.datetime(2010, 1, 1, 11, 0),
         datetime.datetime(2010, 1, 1, 13, 0))

        >>> for x in datetimerange(datetime.datetime(2010, 1, 1, 11),
        ...                        datetime.datetime(2010, 1, 2)):
        ...     x
        ...
        (datetime.datetime(2010, 1, 1, 11, 0), None)

        >>> for x in datetimerange(datetime.datetime(2010, 1, 1, 11),
        ...                        datetime.datetime(2010, 1, 2, 13)):
        ...     x
        ...
        (datetime.datetime(2010, 1, 1, 11, 0), None)
        (datetime.datetime(2010, 1, 2, 0, 0),
         datetime.datetime(2010, 1, 2, 13, 0))

        >>> for x in datetimerange(datetime.datetime(2010, 1, 1, 11),
        ...                        datetime.datetime(2010, 1, 5, 13)):
        ...     x
        ...
        (datetime.datetime(2010, 1, 1, 11, 0), None)
        (datetime.datetime(2010, 1, 2, 0, 0), None)
        (datetime.datetime(2010, 1, 3, 0, 0), None)
        (datetime.datetime(2010, 1, 4, 0, 0), None)
        (datetime.datetime(2010, 1, 5, 0, 0),
         datetime.datetime(2010, 1, 5, 13, 0))

    """
    if start.date() == stop.date():
        yield start, stop
        return
    else:
        yield start, None
        cur = (start.replace(hour=0, minute=0, second=0, microsecond=0) +
               datetime.timedelta(days=1))
        while cur.date() < stop.date():
            yield cur, None
            cur += datetime.timedelta(days=1)
        if cur != stop:
            yield cur, stop
        else:
            return


def _get_or_create_group(file, group):
    """Get or create a group in the datafile"""

    try:
        group = file.get_node(group)
    except tables.NoSuchNodeError:
        parent, newgroup = os.path.split(group)
        file.create_group(parent, newgroup, 'Data group',
                         createparents=True)
    return group


def _get_or_create_node(file, group, src_node):
    """Get or create a node based on a source node"""

    try:
        node = file.get_node(group, src_node.name)
    except tables.NoSuchNodeError:
        if type(src_node) == tables.Table:
            node = file.create_table(group, src_node.name,
                                    src_node.description, src_node.title)
        elif type(src_node) == tables.VLArray:
            node = file.create_vlarray(group, src_node.name, src_node.atom,
                                      src_node.title)
        else:
            raise Exception("Unknown node class: %s" % type(src_node))

    return node
