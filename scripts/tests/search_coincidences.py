#!/usr/bin/env python

"""Search for coincidences

This script tests the process of searching for coincidences.

"""
import datetime

import tables

from sapphire.publicdb import download_data
from sapphire.analysis import coincidences


STATIONS = [501, 503, 506]
START = datetime.datetime(2013, 1, 1)
END = datetime.datetime(2013, 1, 2)


if __name__ == '__main__':
    station_groups = ['/s%d' % u for u in STATIONS]

    data = tables.open_file('data.h5', 'w')
    for station, group in zip(STATIONS, station_groups):
        download_data(data, group, station, START, END)

    coincidences = coincidences.Coincidences(data, '/coincidences',
                                             station_groups)
    coincidences.search_and_store_coincidences()

    # This is the manual method
    #coincidences.search_coincidences()
    #coincidences.process_events(overwrite=True)
    #coincidences.store_coincidences()
