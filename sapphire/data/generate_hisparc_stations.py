""" Generate HiSPARC stations JSON

This script updates the JSON containing a list of HiSPARC stations. For
each station the 'info' object from the API is included. This includes
station numbers, scintillator positions, (sub)cluster name and GPS
locations. If internet is unavailable :mod:`~sapphire.api` uses this
JSON.

Cluster objects track station and detector positions over time. To facilitate
this the gps and station layout data for the HiSPARC Network is stored.

"""
from json import dump
from datetime import date
from os import path, extsep, mkdir

from sapphire.api import Station, Network, LOCAL_BASE, SRC_BASE
from sapphire.utils import pbar


JSON_FILE = path.join(LOCAL_BASE, 'hisparc_stations.json')


def generate_json():
    """Get the API info data for each station"""

    station_numbers = Network().station_numbers()
    station_info = {number: Station(number).info
                    for number in pbar(station_numbers)}

    return station_info


def save_json(data):
    """Overwrite the existing JSON and include the date"""

    data['_info'] = "HiSPARC station info on %s" % date.today()
    with open(JSON_FILE, 'w') as json_file:
        dump(data, json_file, indent=4, sort_keys=True)


def save_tsv():
    """Get location tsv data for all stations"""

    station_numbers = Network().station_numbers()
    for type in ['gps', 'trigger', 'layout']:
        try:
            mkdir(path.join(LOCAL_BASE, type))
        except OSError:
            pass
        for number in pbar(station_numbers):
            url = Station.src_urls[type].format(station_number=number)
            try:
                data = Station._retrieve_url(url, base=SRC_BASE)
            except:
                print 'Failed to get %s data for station %d' % (type, number)
                continue
            data = '\n'.join(d for d in data.split('\n')
                             if len(d) and d[0] != '#')
            if data:
                with open(url.strip('/') + extsep + 'tsv', 'w') as tsvfile:
                    tsvfile.write(data)


if __name__ == '__main__':
    data = generate_json()
    save_json(data)
    save_tsv()
