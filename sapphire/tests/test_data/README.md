How to create the test data
===========================


Notes on recreating publicdb.h5 and publicdb_src.h5
---------------------------------------------------

This test data is simply 5 minutes of real data from station 501 downloaded
from the raw data via the Public Database. Data from 2016-4-21 is used because
it contains many types of data: comparator, config, events, singles, and
weather. Because the download script downloads all blobs for the day, those
need to be removed after downloading, and the file repacked otherwise the file
remains large. This can be done using the following script:

    >>> import tables
    >>> from datetime import datetime
    >>> from sapphire.publicdb import download_data
    >>> with tables.open_file('publicdb_temp.h5', 'w') as data:
    ...     download_data(data, '/station_501', 501, datetime(2016, 4, 21),
    ...                   datetime(2016, 4, 21, 0, 1, 30), get_blobs=True)
    ...     max_blobs_id = data.root.s501.events.col('traces').max()
    ...     data.root.s501.blobs.truncate(max_blobs_id + 1)

    $ ptrepack --complevel 9 --complib blosc publicdb_temp.h5 publicdb_src.h5
    $ rm publicdb_temp.h5
    $ cp publicdb_src.h5 publicdb_src_tmp.h5

    >>> import tables
    >>> from datetime import datetime
    >>> from sapphire.publicdb import _store_data
    >>> f = tables.Filters(complevel=1)
    >>> with tables.open_file('publicdb.h5', 'w', filters=f) as data:
    ...     _store_data(data, '/s501', 'publicdb_src_tmp.h5',
    ...                 datetime(2016, 4, 21), datetime(2016, 4, 21, 0, 1))


Other data
----------

The other data is created by the `esd_load_data.py` script one level up.
