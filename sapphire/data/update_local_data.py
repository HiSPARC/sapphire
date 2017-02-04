""" Update local JSON and TSV data

This script updates the local copies of the JSON and TSV data from the Public
Database API. If internet is unavailable the :mod:`~sapphire.api` uses these
files. The use of local data can also be forced to skip calls to the server or
prevented to require fresh data from the server.

Not all available data is included by default because then the SAPPHiRE package
would become to large. It is possible to add those files after installation.

This can be run by simply running the installed script from the command line::

    $ update_local_data

To make the script show information about what it will do add the help flag::

    $ update_local_data --help

"""
from __future__ import print_function

from json import dump, loads
from os import path, extsep, mkdir, makedirs
from itertools import combinations
import argparse
import warnings

from ..clusters import HiSPARCNetwork
from ..api import API, Network, LOCAL_BASE, SRC_BASE
from ..utils import pbar


def update_local_json(progress=True):
    """Get cluster organisation and basic station JSON data"""

    toplevel_types = ['stations', 'subclusters', 'clusters', 'countries']
    if progress:
        print('Downloading JSONs: %s' % '/'.join(toplevel_types))
    for data_type in pbar(toplevel_types, show=progress):
        update_toplevel_json(data_type)

    for arg_type, data_type in [('stations', 'station_info'),
                                ('subclusters', 'stations_in_subcluster'),
                                ('clusters', 'subclusters_in_cluster'),
                                ('countries', 'clusters_in_country')]:
        if progress:
            print('Downloading JSONs: %s' % data_type)
        update_sublevel_json(arg_type, data_type, progress)


def update_local_tsv(progress=True):
    """Get configuration and calibration TSV data for all stations"""

    station_numbers = Network().station_numbers()

    for data_type in ['gps', 'trigger', 'layout', 'voltage', 'current',
                      'electronics', 'detector_timing_offsets']:
        if progress:
            print('Downloading TSVs: %s' % data_type)
        update_sublevel_tsv(data_type, station_numbers)

    # GPS and layout data should now be up to date, local data can be used
    with warnings.catch_warnings(record=True):
        network = HiSPARCNetwork(force_stale=True)

    for data_type in ['station_timing_offsets']:
        if progress:
            print('Downloading TSVs: %s' % data_type)
        update_subsublevel_tsv(data_type, station_numbers, network)


def update_toplevel_json(data_type):
    url = API.urls[data_type]
    try:
        get_and_store_json(url)
    except Exception:
        print('Failed to get %s data' % data_type)


def update_sublevel_json(arg_type, data_type, progress=True):
    subdir = API.urls[data_type].split('/')[0]
    try:
        mkdir(path.join(LOCAL_BASE, subdir))
    except OSError:
        pass

    url = API.urls[arg_type]
    try:
        numbers = [x['number'] for x in loads(API._retrieve_url(url))]
    except Exception:
        if progress:
            print('Failed to get %s data' % data_type)
        return

    kwarg = API.urls[data_type].split('/')[1].strip('{}')
    for number in pbar(numbers, show=progress):
        url = API.urls[data_type].format(**{kwarg: number, 'year': '',
                                            'month': '', 'day': ''})
        try:
            get_and_store_json(url.strip('/'))
        except Exception:
            if progress:
                print('Failed to get %s data for %s %d' %
                      (data_type, arg_type, number))
            return


def update_sublevel_tsv(data_type, station_numbers, progress=True):
    subdir = API.src_urls[data_type].split('/')[0]
    try:
        mkdir(path.join(LOCAL_BASE, subdir))
    except OSError:
        pass

    for number in pbar(station_numbers, show=progress):
        url = API.src_urls[data_type].format(station_number=number,
                                             year='', month='', day='')
        url = url.strip('/') + '/'
        try:
            get_and_store_tsv(url)
        except Exception:
            if progress and data_type != 'layout':
                print('Failed to get %s for station %d' % (data_type, number))
            continue


def update_subsublevel_tsv(data_type, station_numbers, network, progress=True):
    subdir = API.src_urls[data_type].split('/')[0]
    for number1, number2 in pbar(list(combinations(station_numbers, 2)),
                                 show=progress):
        distance = network.calc_distance_between_stations(number1, number2)
        if distance is None or distance > 1e3:
            continue
        try:
            makedirs(path.join(LOCAL_BASE, subdir, str(number1)))
        except OSError:
            pass
        url = API.src_urls[data_type].format(station_1=number1,
                                             station_2=number2)
        try:
            get_and_store_tsv(url)
        except Exception:
            if progress:
                print('Failed to get %s data for station pair %d-%d' %
                      (data_type, number1, number2))


def get_and_store_json(url):
    data = loads(API._retrieve_url(url))
    json_path = path.join(LOCAL_BASE, url.strip('/') + extsep + 'json')
    with open(json_path, 'w') as jsonfile:
        dump(data, jsonfile, indent=4, sort_keys=True)


def get_and_store_tsv(url):
    data = API._retrieve_url(url, base=SRC_BASE)
    # Strip empty and comment lines
    data = '\n'.join(d for d in data.split('\n') if len(d) and d[0] != '#')
    if data:
        # End with empty newline
        data += '\n'
        tsv_path = path.join(LOCAL_BASE, url.strip('/') + extsep + 'tsv')
        with open(tsv_path, 'w') as tsvfile:
            tsvfile.write(data)


def main():
    descr = """Update commonly used local data. This allows the usage of
               api and cluster objects without an internet connection. Or with
               internet the usage of local data can be forced to speed up the
               retrieval of data. This data is already included in SAPPHiRE,
               but this script makes the data up to date."""
    parser = argparse.ArgumentParser(description=descr)
    parser.parse_args()
    update_local_json()
    update_local_tsv()
