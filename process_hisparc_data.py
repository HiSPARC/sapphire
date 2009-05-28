""" Store data from the HiSPARC eventwarehouse using pytables

    This module processes data read from the event and eventdata tables in
    the eventwarehouse database. Timestamps are calculated from the date
    and time columns (remember, in GPS time!). The eventdata columns are
    processed and stored alongside the event data.

"""
import datetime
import time
import os

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class IntegrityError(Error):
    """Exception raised for data integrity errors.

    Attributes:
        message --- error message

    """
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


def process_hisparc_events(events, eventdata, table):
    """Do the actual data processing and storing

    This function concurrently reads the events and eventdata lists and
    builds up the event rows. When a row is complete, i.e. there are no
    more eventdata values for the current event row, it is stored in a
    pytables table.

    This is one nice long algorithmic hack. This is the price we have to
    pay to gain a speed increase of a factor of 20 or so.

    Arguments:
    events          contents from the eventwarehouse event table
    eventdata       contents from the eventwarehouse eventdata table
    table           the destination table

    """
    tablerow = table.row
    data_idx = 0

    # First, make sure we have no 'old' eventdata records in the first few
    # rows.
    event_id = events[0][0]
    try:
        while eventdata[data_idx][0] != event_id:
            data_idx += 1
    except IndexError:
        # We've exhausted all eventdata records, there is no match!
        raise IntegrityError("Eventdata records don't match event " \
                             "records.")

    for row in events:
        # We process the events row by row
        event_id = row[0]
        date = row[1]
        timedelta = row[2]
        nanoseconds = row[3]

        # calculate the timestamp (in GPS time, not UTC time!)
        t = datetime.datetime.combine(date, datetime.time()) + timedelta
        timestamp = time.mktime(t.utctimetuple())

        tablerow['event_id'] = event_id
        tablerow['timestamp'] = timestamp
        tablerow['nanoseconds'] = nanoseconds
        tablerow['ext_timestamp'] = timestamp * 1e9 + nanoseconds

        # get default values for the eventdata
        data = {}
        data['pulseheights'] = tablerow['pulseheights']
        data['integrals'] = tablerow['integrals']

        while True:
            # now process the eventdata row by row, using the current index
            # data_idx
            try:
                data_row = eventdata[data_idx]
            except IndexError:
                # We've exhausted all eventdata. Break to store the current
                # event. Hopefully, we've exhausted events as well.
                break
            if data_row[0] == event_id:
                # the eventdata matches the current event, make sure to
                # read the next eventdata row next time and process the
                # current row
                data_idx += 1
                uploadcode = data_row[1]
                value = data_row[2]

                if uploadcode[:2] == 'PH':
                    key = 'pulseheights'
                elif uploadcode[:2] == 'IN':
                    key = 'integrals'
                else:
                    continue
                idx = int(uploadcode[2]) - 1 
                data[key][idx] = value
            else:
                # The eventdata is not matching this event. Probably we've
                # exhausted this events' eventdata. Break to store this
                # event.
                break
        if (data['pulseheights'] == tablerow['pulseheights']).all() or \
           (data['integrals'] == tablerow['integrals']).all():
            raise IntegrityError("Possibly missing event records.")

        tablerow['pulseheights'] = data['pulseheights']
        tablerow['integrals'] = data['integrals']
        tablerow.append()
        # continue on to the next event

    table.flush()

if __name__ == '__main__':
    print 'Setting time zone to UTC'
    os.environ['TZ'] = 'UTC'
    time.tzset()
