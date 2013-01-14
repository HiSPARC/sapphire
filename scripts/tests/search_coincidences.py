#!/usr/bin/env python

"""Search for coincidences

This script tests the process of searching for coincidences.

"""
import datetime

import tables

from sapphire.publicdb import download_data
from sapphire.analysis import coincidences


STATIONS = [501, 503, 506]


if __name__ == '__main__':
    station_groups = ['/s%d' % u for u in STATIONS]

    data = tables.openFile('testdata.h5', 'a')
    for station, group in zip(STATIONS, station_groups):
        if group not in data:
            download_data(data, group, station,
                          datetime.datetime(2013, 1, 1),
                          datetime.datetime(2013, 1, 2), get_blobs=False)

    coincidences = coincidences.Coincidences(data, '/coincidences',
                                             station_groups,
                                             overwrite=True)
    coincidences.search_and_store_coincidences()

    # This is the manual method
    #coincidences.search_coincidences()
    #coincidences.process_events(overwrite=True)
    #coincidences.store_coincidences()
