""" Add more local JSON and TSV data

Add additional local data, to be used by :mod:`~sapphire.api` if internet is
unavailable. The use of local data can also be forced to skip calls to the
server or prevented to require fresh data from the server.

This data is not included by default because then the SAPPHiRE package would
become to large. By running this script the data is added after installation.

"""
from os import path, extsep, mkdir
import argparse

from ..api import API, Network, LOCAL_BASE, SRC_BASE
from ..utils import pbar


def update_additional_local_tsv():
    """Get location tsv data for all stations"""

    station_numbers = Network().station_numbers()

    for type in ['eventtime']:
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


def main():
    descr = """Add additional data to local data, or update already downloaded
             data. Making data available locally can greatly speed up a
             program which uses this data. Approximately 100 MB of data will
             be downloaded. The data contains the eventtime data, i.e. hourly
             number of events for all stations."""
    parser = argparse.ArgumentParser(description=descr)
    parser.parse_args()
    update_additional_local_tsv()
