""" Update local JSON and TSV data

This script updates the local copies of the JSON and TSV data from the Public
Database API. If internet is unavailable the :mod:`~sapphire.api` uses these
files. The use of local data can also be forced to skip calls to the server or
prevented to require fresh data from the server.

Not all available data is included by default because then the SAPPHiRE package
would become to large. It is possible to add those files after installation.

"""
from json import dump, loads
from os import path, extsep, mkdir

from sapphire.api import API, Network, LOCAL_BASE, SRC_BASE
from sapphire.utils import pbar


def update_local_json():
    for type in pbar(['stations', 'subclusters', 'clusters', 'countries']):
        url = API.urls[type]
        try:
            data = loads(API._retrieve_url(url))
        except:
            print 'Failed to get %s data' % type
            continue
        if data:
            json_path = path.join(LOCAL_BASE,
                                  url.strip('/') + extsep + 'json')
            with open(json_path, 'w') as jsonfile:
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
        for number in pbar(numbers):
            url = API.urls[type].format(**{kwarg: number,
                                           'year': '', 'month': '', 'day': ''})
            try:
                data = loads(API._retrieve_url(url.strip('/')))
            except:
                print 'Failed to get %s data for %s %d' % (type, arg, number)
                continue
            if data:
                json_path = path.join(LOCAL_BASE,
                                      url.strip('/') + extsep + 'json')
                with open(json_path, 'w') as jsonfile:
                    dump(data, jsonfile, indent=4, sort_keys=True)


def update_local_tsv():
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
                if type != 'layout':
                    print 'Failed to get %s for station %d' % (type, number)
                continue
            data = '\n'.join(d for d in data.split('\n')
                             if len(d) and d[0] != '#')
            if data:
                tsv_path = path.join(LOCAL_BASE,
                                     url.strip('/') + extsep + 'tsv')
                with open(tsv_path, 'w') as tsvfile:
                    tsvfile.write(data)


if __name__ == '__main__':
    update_local_json()
    update_local_tsv()
