#!/usr/bin/env python

"""Process HiSPARC events without traces

This script tests the ProcessEventsWithoutTraces class, as well as the
ProcessIndexedEventsWithoutTraces class.

"""
import datetime

import tables

from sapphire.publicdb import download_data
from sapphire.analysis import process_events


if __name__ == '__main__':
    data = tables.open_file('testdata.h5', 'a')
    if '/s501' not in data:
        download_data(data, '/s501', 501, datetime.datetime(2013, 1, 1),
                      datetime.datetime(2013, 1, 2), get_blobs=False)
    process = process_events.ProcessEventsWithoutTraces(data, '/s501')
    process.process_and_store_results(overwrite=True)
    offsets = process.determine_detector_timing_offsets()
    print "Offsets:", offsets
