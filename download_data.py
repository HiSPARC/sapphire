import sys
import os.path
import tables

from create_tables import create_tables
from get_hisparc_data import get_hisparc_data
from process_hisparc_data import process_hisparc_events

def download_data(datafile, chunksize=50000, limit=None):
    """ Download HiSPARC data and append it to the data file

    This function downloads HiSPARC data, continuing with new data and
    appending it to the pytables data file.

    Arguments:
    datafile            the pytables data file
    chunksize           number of events at a time which are downloaded and
                        processed

    """
    table = datafile.root.hisparc.events

    offset = len(table)
    count = 0
    while True:
        print "Downloading %d events, starting from offset %d... " % \
            (chunksize, offset)
        events, eventdata = get_hisparc_data(limit=chunksize,
                                             offset=offset)
        print "done."

        if events:
            print "Processing events... ",
            sys.stdout.flush()
            process_hisparc_events(events, eventdata, table)
            offset += chunksize
            print "done."
        else:
            break

        count += 1
        if limit and count >= limit:
            print "Limit reached, breaking."
            break

if __name__ == '__main__':
    if not os.path.isfile('data_new.h5'):
        create_tables()

    datafile = tables.openFile('data_new.h5', 'a')
    download_data(datafile)
