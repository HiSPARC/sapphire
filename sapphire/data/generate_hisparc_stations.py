from json import dump
from datetime import date
from os import path

from sapphire.api import Network, Station
from sapphire.utils import pbar


JSON_FILE = path.join(path.dirname(__file__), 'hisparc_stations.json')


def generate_json():
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
    data['_info'] = "HiSPARC station info on %s" % date.today()
    with open(JSON_FILE, 'w') as json_file:
        dump(data, json_file, indent=4, sort_keys=True)


if __name__ == '__main__':
    data = generate_json()
    save_json(data)
