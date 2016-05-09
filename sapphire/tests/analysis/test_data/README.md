How to create the test data
===========================


Notes on recreating process_events.h5
----------------------------------

This test data is simply 5 minutes of real data from station 501
downloaded from the raw data via the Public Database.
Because the download script downloads all blobs those need to be removed
after downloading, and the file repacked otherwise you are left with
60MB of data. This can be done using the following script:

    >>> import tables
    >>> from datetime import datetime
    >>> from sapphire.publicdb import download_data
    >>> with tables.open_file('process_events_temp.h5', 'w') as data:
    ...     download_data(data, '/s501', 501, datetime(2010, 9, 1),
    ...                   datetime(2010, 9, 1, 0, 5), get_blobs=True)
    ...     max_trace_id = data.root.s501.events.col('traces').max()
    ...     data.root.s501.blobs.truncate(max_trace_id + 1)

    $ ptrepack --complevel 9 --complib blosc process_events_temp.h5 process_events.h5
    $ rm process_events_temp.h5


Notes on recreating coincidences.h5
---------------------------------------

    >>> from mock import patch
    >>> import tables
    >>> from sapphire import download_data, Coincidences
    >>> from datetime import datetime
    >>> filters = tables.Filters(complevel=1)
    >>> start = datetime(2012, 1, 1, 0, 0, 0)
    >>> end = datetime(2012, 1, 1, 0, 2, 0)
    >>> with tables.open_file('coincidences.h5', 'w', filters=filters) as data:
    ...     download_data(data, '/station_501', 501, start, end, progress=False)
    ...     download_data(data, '/station_502', 502, start, end, progress=False)
    ...     with patch('sapphire.analysis.process_events.ProcessIndexedEventsWithoutTraces'):
    ...         c = Coincidences(data, '/coincidences', ['/station_501', '/station_502'], progress=False)
    ...         c.search_and_store_coincidences()


Notes on recreating esd_coincidences.h5
---------------------------------------

    >>> import tables
    >>> from sapphire import download_data, CoincidencesESD
    >>> from datetime import datetime
    >>> filters = tables.Filters(complevel=1)
    >>> start = datetime(2012, 1, 1, 0, 0, 0)
    >>> end = datetime(2012, 1, 1, 0, 2, 0)
    >>> with tables.open_file('esd_coincidences.h5', 'w', filters=filters) as data:
    ...     download_data(data, '/station_501', 501, start, end, progress=False)
    ...     download_data(data, '/station_502', 502, start, end, progress=False)
    ...     c = CoincidencesESD(data, '/coincidences', ['/station_501', '/station_502'], progress=False)
    ...     c.search_and_store_coincidences(station_numbers=[501, 502])
