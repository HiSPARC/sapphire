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
from multiprocessing import Process, Queue, Event
import signal

from hisparc.eventwarehouse import get_events, process_events

interrupt = Event()


def handler(signum, frame):
    """Handler to signal a KeyboardInterrupt to all threads"""
    if not interrupt.is_set():
        print "KeyboardInterrupt, please wait..."
        interrupt.set()
    else:
        pass

def download(queue, station_id, chunksize, offset, limit, get_traces):
    """Download HiSPARC data

    This function downloads HiSPARC data, starting from an offset.  It
    downloads chunksize events at a time and ends after a limit number of
    downloads.  Downloaded data will be put in a Queue object.

    :param queue: a multiprocessing Queue object
    :param chunksize: number of events at a time which are downloaded and
        processed
    :param offset: the number of events which should be skipped
    :param limit: the maximum number of chunks to process
    :param get_traces: boolean, select whether traces should be fetched

    """
    count = 0
    # Loop until interrupt is signalled
    while not interrupt.is_set():
        print "Downloading %d events, starting from offset %d... " % \
            (chunksize, offset)
        events, eventdata, calculateddata = get_events(station_id,
                                                       limit=chunksize,
                                                       offset=offset,
                                                       get_traces=get_traces)
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

    if interrupt.is_set():
        print ("Downloader interrupted, please wait for processor to "
               "finish...")

def process(queue, eventstable, traces):
    """Append HiSPARC data to the data file

    This function processes the data which has been put in Queue object,
    appending it to the pytables data table.

    :param queue: a multiprocessing Queue object
    :param events: the PyTables destination events table
    :param traces: the PyTables destination traces array

    """

    while True:
        try:
            events, eventdata, calculateddata = queue.get()
        except IOError:
            # Probably due to repeatedly pressing Ctrl-C, restart call
            continue

        if events:
            print "Processing events... "
            process_events(events, eventdata, calculateddata,
                           eventstable, traces)
            print "done."
        else:
            print "No more events, shutting down."
            break


def start_download(file, group, station_id=601, limit=1, chunksize=5000,
                   get_traces=False):
    """Start a multi-process download

    A convenience function to start a download and process the data in
    separate processes to speed up download on multi-core machines.

    :param file: The PyTables datafile handler
    :param group: The PyTables destination group, which should have an
        events table and traces array
    :param station_id: The HiSPARC station number for which to get events
    :param limit: The number of chunks to download
    :param chunksize: The number of events in one chunk
    :param get_traces: boolean, select whether traces should be fetched

    Example usage::

        >>> import tables
        >>> data = tables.openFile('test.h5', 'a')
        >>> from hisparc.utils import download_data
        >>> download_data.start_download(data,
        ... '/hisparc/sciencepark/station501', limit=2, get_traces=True)
        Starting download subprocess...
        Waiting for download to shut down...
        Downloading 5000 events, starting from offset 0... 
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
    # create a custom signal handler to elegantly handle KeyboardInterrupt
    # in all handlers
    signal.signal(signal.SIGINT, handler)
    signal.siginterrupt(signal.SIGINT, False)
    interrupt.clear()

    events = file.getNode(group, 'events')
    traces = file.getNode(group, 'traces')
    offset = len(events)

    queue = Queue(maxsize=2)
    downloader = Process(target=download, args=(queue, station_id),
                         kwargs={'offset': offset, 'chunksize': chunksize,
                                 'limit': limit, 'get_traces': get_traces})

    print "Starting download subprocess..."
    downloader.start()

    print "Waiting for download to shut down..."
    process(queue, events, traces)
    downloader.join()
