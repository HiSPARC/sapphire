""" Easy to use HiSPARC data downloader and processor

    This script downloads and processes HiSPARC event data.  It uses the
    multiprocessing module to split up the process in two subprocesses, a
    downloader and a processor, using a Queue object to manage the data.
    PyTables is used to store the data in a HDF5 file.  Not much of the
    functionality is actually in this script, but is rather imported from
    other scripts.  Therefore, this script is just a convenience.

"""
import os.path
import tables
from multiprocessing import Process, Queue

from create_tables import create_tables
from get_hisparc_data import get_hisparc_data
from process_hisparc_events import process_hisparc_events

def downloader(queue, chunksize=50000, offset=0, limit=None):
    """Download HiSPARC data

    This function downloads HiSPARC data, starting from an offset.  It
    downloads chunksize events at a time and ends after a limit number of
    downloads.  Downloaded data will be put in a Queue object.

    Arguments:
    queue               a multiprocessing Queue object
    chunksize           number of events at a time which are downloaded and
                        processed
    offset              the number of events which should be skipped
    limit               the maximum number of chunks to process

    """
    count = 0
    try:
        while True:
            print "Downloading %d events, starting from offset %d... " % \
                (chunksize, offset)
            events, eventdata = get_hisparc_data(limit=chunksize,
                                                 offset=offset)
            print "done."

            if events:
                print "Putting events in queue..."
                queue.put((events, eventdata))
                offset += chunksize
            else:
                print "No more events, shutting down."
                break

            count += 1
            if limit and count >= limit:
                print "Limit reached, shutting down."
                break

        # Signalling shut down
        queue.put((None, None))
    except KeyboardInterrupt:
        print "Downloader received KeyboardInterrupt, shutting down..."

def processor(queue, table):
    """Download HiSPARC data and append it to the data file

    This function processes the data which has been put in Queue object,
    appending it to the pytables data table.

    Arguments:
    queue               a multiprocessing Queue object
    table               the pytables data table

    """
    table = datafile.root.hisparc.events

    try:
        while True:
            events, eventdata = queue.get()

            if events:
                print "Processing events... "
                process_hisparc_events(events, eventdata, table)
                print "done."
            else:
                print "No more events, shutting down."
                break
    except KeyboardInterrupt:
        print "Processor received KeyboardInterrupt, shutting down..."


if __name__ == '__main__':
    if not os.path.isfile('data_new.h5'):
        create_tables()

    datafile = tables.openFile('data_new.h5', 'a')
    table =  datafile.root.hisparc.events
    offset = len(table)

    queue = Queue(maxsize=2)
    downloader = Process(target=downloader, args=(queue,),
                         kwargs={'offset': offset})
    processor = Process(target=processor, args=(queue, table))

    print "Starting subprocesses..."
    downloader.start()
    processor.start()

    print "Waiting for subprocesses to shut down..."
    try:
        downloader.join()
        processor.join()
    except KeyboardInterrupt:
        print "Main program received KeyboardInterrupt."
