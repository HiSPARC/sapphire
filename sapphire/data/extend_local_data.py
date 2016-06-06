""" Add more local JSON and TSV data

Add additional local data, to be used by :mod:`~sapphire.api` if internet is
unavailable. The use of local data can also be forced to skip calls to the
server or prevented to require fresh data from the server.

This data is not included by default because then the SAPPHiRE package would
become to large. By running this script the data is added after installation.

"""
from os import path, extsep, mkdir
from itertools import combinations

from sapphire import HiSPARCNetwork
from sapphire.api import API, Network, LOCAL_BASE, SRC_BASE
from sapphire.utils import pbar


def update_additional_local_tsv():
    """Get location tsv data for all stations"""

    station_numbers = Network().station_numbers()

    for type in ['eventtime', 'detector_timing_offsets']:
        try:
            mkdir(path.join(LOCAL_BASE, type))
        except OSError:
            pass
        for number in pbar(station_numbers):
            url = API.src_urls[type].format(station_number=number,
                                            year='', month='', day='', hour='')
            try:
                data = API._retrieve_url(url.strip('/'), base=SRC_BASE)
            except:
                print 'Failed to get %s data for station %d' % (type, number)
                continue
            data = '\n'.join(d for d in data.split('\n')
                             if len(d) and d[0] != '#')
            if data:
                tsv_path = path.join(LOCAL_BASE,
                                     url.strip('/') + extsep + 'tsv')
                with open(tsv_path, 'w') as tsvfile:
                    tsvfile.write(data)

    type = 'station_timing_offsets'
    network = HiSPARCNetwork()

    try:
        mkdir(path.join(LOCAL_BASE, type))
    except OSError:
        pass
    for number1, number2 in pbar(combinations(station_numbers, 2)):
        distance = network.calc_distance_between_stations(number1, number2)
        if distance is None or distance > 1e3:
            continue
        try:
            mkdir(path.join(LOCAL_BASE, type, str(number1)))
        except OSError:
            pass
        url = API.src_urls[type].format(station_1=number1, station_2=number2)
        try:
            data = API._retrieve_url(url.strip('/'), base=SRC_BASE)
        except:
            print ('Failed to get %s data for station pair %d-%d' %
                   (type, number1, number2))
            continue
        data = '\n'.join(d for d in data.split('\n') if len(d) and d[0] != '#')
        if data:
            tsv_path = path.join(LOCAL_BASE, url.strip('/') + extsep + 'tsv')
            with open(tsv_path, 'w') as tsvfile:
                tsvfile.write(data)


if __name__ == '__main__':
    update_additional_local_tsv()
