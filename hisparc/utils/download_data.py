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

from hisparc.eventwarehouse import get_events, process_events

def download(queue, chunksize=50000, offset=0, limit=None):
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
            events, eventdata, calculateddata = get_events(601,
                                                           limit=chunksize,
                                                           offset=offset)
            print "done."

            if events:
                print "Putting events in queue..."
                queue.put((events, eventdata, calculateddata))
                offset += chunksize
            else:
                print "No more events, shutting down."
                break

            count += 1
            if limit and count >= limit:
                print "Limit reached, shutting down."
                break

        # Signalling shut down
        queue.put((None, None, None))
    except KeyboardInterrupt:
        print "Downloader received KeyboardInterrupt, shutting down..."

def process(queue, datafile):
    """Append HiSPARC data to the data file

    This function processes the data which has been put in Queue object,
    appending it to the pytables data table.

    Arguments:
    queue               a multiprocessing Queue object
    datafile            the pytables data file

    """
    try:
        while True:
            events, eventdata, calculateddata = queue.get()

            if events:
                print "Processing events... "
                process_events(events, eventdata, calculateddata,
                               datafile.root.hisparc.events,
                               datafile.root.hisparc.traces)
                print "done."
            else:
                print "No more events, shutting down."
                break
    except KeyboardInterrupt:
        print "Processor received KeyboardInterrupt, shutting down..."


def start_download(filename, limit=1, chunksize=5000):
    """Start a multi-process download

    A convenience function to start a download and process the data in
    separate processes to speed up download on multi-core machines.

    Arguments:
    filename            The filename of the datafile
    limit               The number of chunks to download
    chunksize           The number of events in one chunk

    """
    if not os.path.isfile(filename):
        create_tables(filename)

    datafile = tables.openFile(filename, 'a')
    offset = len(datafile.root.hisparc.events)

    queue = Queue(maxsize=2)
    downloader = Process(target=download, args=(queue,),
                         kwargs={'offset': offset, 'chunksize': chunksize,
                                 'limit': limit})
    processor = Process(target=process, args=(queue, datafile))

    print "Starting subprocesses..."
    downloader.start()
    processor.start()

    print "Waiting for subprocesses to shut down..."
    try:
        downloader.join()
        processor.join()
    except KeyboardInterrupt:
        print "Main program received KeyboardInterrupt."
        downloader.join()
        processor.join()

    datafile.close()
