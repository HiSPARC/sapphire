import sys
import os.path
import calendar
import base64
import tables
import logging
import datetime

import storage
from upload_codes import eventtype_upload_codes

logger = logging.getLogger('writer.store_events')


def store_event(datafile, cluster, station_id, event):
    """Stores an event in the h5 filesystem

    :param datafile: the h5 data file
    :param cluster: the name of the cluster to which the station belongs
    :param station_id: the id of the station this event belongs to
    :param event: the event to store

    """
    eventheader = event['header']
    eventdatalist = event['datalist']
    eventtype = eventheader['eventtype_uploadcode']

    try:
        upload_codes = eventtype_upload_codes[eventtype]
    except KeyError:
        logger.error('Unknown event type: %s, discarding event' %
                     eventtype)
        return

    parentnode = storage.get_or_create_station_group(datafile, cluster,
                                                     station_id)
    table = storage.get_or_create_node(datafile, parentnode,
                                       upload_codes['_tablename'])
    blobs = storage.get_or_create_node(datafile, parentnode, 'blobs')

    row = table.row
    row['event_id'] = table.nrows + 1
    # make a unix-like timestamp
    timestamp = calendar.timegm(eventheader['datetime'].utctimetuple())
    nanoseconds = eventheader['nanoseconds']
    # make an extended timestamp, which is the number of nanoseconds since
    # epoch
    ext_timestamp = timestamp * int(1e9) + nanoseconds
    row['timestamp'] = timestamp

    if upload_codes['_has_ext_time']:
        # This is e.g. a HiSPARC coincidence or comparator message,
        # extended timing information is available
        row['nanoseconds'] = nanoseconds
        row['ext_timestamp'] = ext_timestamp

    # get default values for the data
    data = {}
    for key, value in upload_codes.items():
        if key[0] != '_':
            # private meta information starts with a _ (e.g. _tablename)
            data[key] = row[value]

    # process event data
    for item in eventdatalist:
        # uploadcode: EVENTRATE, PH1, IN3, etc.
        uploadcode = item['data_uploadcode']
        # value: actual data value
        value = item['data']

        if data_is_blob(uploadcode, upload_codes['_blobs']):
            # data should be stored inside the blob array, ...
            if uploadcode[:-1] == 'TR':
                # traces are base64 encoded
                value = base64.decodestring(value)
            blobs.append(value)
            # ... with a pointer stored in the event table
            value = len(blobs) - 1

        if uploadcode[-1].isdigit():
            # uploadcode: PH1, IN3, etc.
            key, index = uploadcode[:-1], int(uploadcode[-1]) - 1
            if key in data:
                data[key][index] = value
            else:
                logger.warning('Datatype not known on server side: %s '
                               '(%s)' % (key, eventtype))
        else:
            # uploadcode: EVENTRATE, RED, etc.
            if uploadcode in data:
                data[uploadcode] = value
            else:
                logger.warning('Datatype not known on server side: %s '
                               '(%s)' % (uploadcode, eventtype))

    # write data values to row
    for key, value in upload_codes.items():
        if key[0] != '_':
            # private meta information starts with a _ (e.g. _tablename)
            row[value] = data[key]

    row.append()
    table.flush()
    blobs.flush()


def data_is_blob(uploadcode, blob_types):
    """Determine if data is a variable length binary value (blob)"""

    if uploadcode[-1].isdigit():
        if uploadcode[:-1] in blob_types:
            return True
    elif uploadcode in blob_types:
        return True
    return False


def store_event_list(data_dir, station_id, cluster, event_list):
    """Store a list of events"""

    prev_date = None
    datafile = None
    for event in event_list:
        date = event['header']['datetime'].date()
        if date != prev_date:
            if datafile:
                datafile.close()
            datafile = storage.open_or_create_file(data_dir, date)
            prev_date = date

        store_event(datafile, cluster, station_id, event)

    if datafile:
        datafile.close()
