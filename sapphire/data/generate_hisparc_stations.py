""" Generate HiSPARC stations JSON

This script updates the JSON containing a list of HiSPARC stations. For
each station the 'info' object from the API is included. This includes
station numbers, scintillator positions, (sub)cluster name and GPS
locations. If internet is unavailable :mod:`~sapphire.api` uses this
JSON.

Cluster objects track station and detector positions over time. To facilitate
this the gps and station layout data for the HiSPARC Network is stored.

"""
from json import dump, loads
from datetime import date
from os import path, extsep, mkdir

from sapphire.api import API, Station, Network, LOCAL_BASE, API_BASE, SRC_BASE
from sapphire.utils import pbar


def save_json():
    for type in ['stations', 'subclusters', 'clusters', 'countries']:
        url = API.urls[type]
        try:
            data = loads(API._retrieve_url(url))
        except:
            print 'Failed to get %s data' % type
            continue
        if data:
            with open(url.strip('/') + extsep + 'json', 'w') as jsonfile:
                dump(data, jsonfile, indent=4, sort_keys=True)

    for arg, kwarg, type in [
            ('stations', 'station_number', 'station_info'),
            ('subclusters', 'subcluster_number', 'stations_in_subcluster'),
            ('clusters', 'cluster_number', 'subclusters_in_cluster'),
            ('countries', 'country_number', 'clusters_in_country')]:
        try:
            if arg == 'stations':
                mkdir(path.join(LOCAL_BASE, arg[:-1]))
            else:
                mkdir(path.join(LOCAL_BASE, arg))
        except OSError:
            pass
        url = API.urls[arg]
        try:
            numbers = [x['number'] for x in loads(API._retrieve_url(url))]
        except:
            print 'Failed to get %s data' % type
            continue
        for number in numbers:
            url = API.urls[type].format(**{kwarg: number,
                                           'year': '', 'month': '', 'day': ''})
            try:
                data = loads(API._retrieve_url(url.strip('/')))
            except:
                print 'Failed to get %s data for %s %d' % (type, arg, number)
                continue
            if data:
                with open(url.strip('/') + extsep + 'json', 'w') as jsonfile:
                    dump(data, jsonfile, indent=4, sort_keys=True)


def save_tsv():
    """Get location tsv data for all stations"""

    station_numbers = Network().station_numbers()
    for type in ['gps', 'trigger', 'layout', 'voltage', 'current']:
        try:
            mkdir(path.join(LOCAL_BASE, type))
        except OSError:
            pass
        for number in pbar(station_numbers):
            url = API.src_urls[type].format(station_number=number)
            try:
                data = API._retrieve_url(url, base=SRC_BASE)
            except:
                print 'Failed to get %s data for station %d' % (type, number)
                continue
            data = '\n'.join(d for d in data.split('\n')
                             if len(d) and d[0] != '#')
            if data:
                with open(url.strip('/') + extsep + 'tsv', 'w') as tsvfile:
                    tsvfile.write(data)


if __name__ == '__main__':
    save_json()
    save_tsv()
