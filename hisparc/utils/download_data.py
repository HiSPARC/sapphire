""" Easy to use HiSPARC data downloader and processor

    This script downloads and processes HiSPARC event data.  It uses the
    multiprocessing module to split up the process in two subprocesses, a
    downloader and a processor, using a Queue object to manage the data.
    PyTables is used to store the data in a HDF5 file.  Not much of the
    functionality is actually in this script, but is rather imported from
    other scripts.  Therefore, this script is just a convenience.

    You want to use the :func:`start_download` function.

"""
import os.path
import tables
from multiprocessing import Process, Queue

from create_tables import create_tables

from hisparc.eventwarehouse import get_events, process_events

def download(queue, station_id, chunksize=50000, offset=0, limit=None):
    """Download HiSPARC data

    This function downloads HiSPARC data, starting from an offset.  It
    downloads chunksize events at a time and ends after a limit number of
    downloads.  Downloaded data will be put in a Queue object.

    :param queue: a multiprocessing Queue object
    :param chunksize: number of events at a time which are downloaded and
        processed
    :param offset: the number of events which should be skipped
    :param limit: the maximum number of chunks to process

    """
    count = 0
    try:
        while True:
            print "Downloading %d events, starting from offset %d... " % \
                (chunksize, offset)
            events, eventdata, calculateddata = get_events(station_id,
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

    :param queue: a multiprocessing Queue object
    :param datafile: the pytables data file

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


def start_download(filename, station_id=601, limit=1, chunksize=5000):
    """Start a multi-process download

    A convenience function to start a download and process the data in
    separate processes to speed up download on multi-core machines.

    :param filename: The filename of the datafile
    :param limit: The number of chunks to download
    :param chunksize: The number of events in one chunk

    Example usage::

        >>> from hisparc.utils import download_data
        >>> download_data.start_download('testdata.h5', limit=2)
        Creating a new PyTables data file...  done.
        Starting subprocesses...
        Downloading 5000 events, starting from offset 0... 
        Waiting for subprocesses to shut down...
        Time window:  2008-07-01 14:29:45 2008-07-01 14:52:44
        done.
        Putting events in queue...
        Downloading 5000 events, starting from offset 5000... 
        Processing events... 
        Time window:  2008-07-01 14:52:44 2008-07-01 15:15:18
        done.
        Putting events in queue...
        Limit reached, shutting down.
        done.
        Processing events... 
        done.
        No more events, shutting down.

    """
    if not os.path.isfile(filename):
        create_tables(filename)

    datafile = tables.openFile(filename, 'a')
    offset = len(datafile.root.hisparc.events)

    queue = Queue(maxsize=2)
    downloader = Process(target=download, args=(queue, station_id),
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
