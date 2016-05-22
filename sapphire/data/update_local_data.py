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
        update_toplevel_json(type)

    for arg_type, type in [('stations', 'station_info'),
                           ('subclusters', 'stations_in_subcluster'),
                           ('clusters', 'subclusters_in_cluster'),
                           ('countries', 'clusters_in_country')]:
        update_sublevel_json(arg_type, type)


def update_toplevel_json(type):
    url = API.urls[type]
    try:
        get_and_store_json(url)
    except:
        print 'Failed to get %s data' % type


def update_sublevel_json(arg_type, type):
    subdir = API.urls[type].split('/')[0]
    try:
        mkdir(path.join(LOCAL_BASE, subdir))
    except OSError:
        pass

    url = API.urls[arg_type]
    try:
        numbers = [x['number'] for x in loads(API._retrieve_url(url))]
    except:
        print 'Failed to get %s data' % type
        return

    kwarg = API.urls[type].split('/')[1].strip('{}')
    for number in pbar(numbers):
        url = API.urls[type].format(**{kwarg: number,
                                       'year': '', 'month': '', 'day': ''})
        try:
            get_and_store_json(url.strip('/'))
        except:
            print 'Failed to get %s data for %s %d' % (type, arg_type, number)
            return


def update_local_tsv():
    """Get configuration tsv data for all stations"""

    station_numbers = Network().station_numbers()

    for type in ['gps', 'trigger', 'layout', 'voltage', 'current',
                 'electronics']:
        subdir = API.src_urls[type].split('/')[0]
        try:
            mkdir(path.join(LOCAL_BASE, subdir))
        except OSError:
            pass

        for number in pbar(station_numbers):
            url = API.src_urls[type].format(station_number=number)
            try:
                get_and_store_tsv(url)
            except:
                if type != 'layout':
                    print 'Failed to get %s for station %d' % (type, number)
                continue


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
        tsv_path = path.join(LOCAL_BASE, url.strip('/') + extsep + 'tsv')
        with open(tsv_path, 'w') as tsvfile:
            tsvfile.write(data)


if __name__ == '__main__':
    update_local_json()
    update_local_tsv()
