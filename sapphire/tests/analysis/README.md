How to create the test data
===========================


Notes on recreating PE-testdata.h5
----------------------------------

This test data is simply 5 minutes of real data from station 501
downloaded from the raw data via the Public Database.
Because the download script downloads all blobs those need to be removed
after downloading, and the file repacked otherwise you are left with
60MB of data. This can be done using the following script:

    >>> import tables
    >>> from datetime import datetime
    >>> from sapphire.publicdb import download_data
    >>> with tables.openFile('PE-testdata_temp.h5', 'w') as data:
    ...     download_data(data, '/s501', 501, datetime(2010, 9, 1),
    ...                   datetime(2010, 9, 1, 0, 5), get_blobs=True)
    ...     max_trace_id = data.root.s501.events.col('traces').max()
    ...     data.root.s501.blobs.truncate(max_trace_id + 1)

    $ ptrepack --complevel 9 --complib blosc PE-testdata_temp.h5 PE-testdata.h5
    $ rm PE-testdata_temp.h5
