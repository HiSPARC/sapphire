""" Generate HiSPARC stations JSON

This script updates the JSON containing a list of HiSPARC stations. For
each station the 'info' object from the API is included. This includes
station numbers, scintillator positions, (sub)cluster name and GPS
locations. If internet is unavailable :mod:`~sapphire.api` uses this
JSON.

"""
from json import dump
from datetime import date
from os import path

from sapphire.api import Network, Station
from sapphire.utils import pbar


JSON_FILE = path.join(path.dirname(__file__), 'hisparc_stations.json')


def generate_json():
    """Get the API info data for each station"""

    station_numbers = Network().station_numbers()
    station_info = {}

    for number in pbar(station_numbers):
        try:
            station = Station(number)
            station_info[number] = station.info
        except:
            continue

    return station_info


def save_json(data):
    """Overwrite the existing JSON and include the date"""

    data['_info'] = "HiSPARC station info on %s" % date.today()
    with open(JSON_FILE, 'w') as json_file:
        dump(data, json_file, indent=4, sort_keys=True)


if __name__ == '__main__':
    data = generate_json()
    save_json(data)
