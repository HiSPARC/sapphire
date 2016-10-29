""" Add more local JSON and TSV data

Add additional local data, to be used by :mod:`~sapphire.api` if internet is
unavailable. The use of local data can also be forced to skip calls to the
server or prevented to require fresh data from the server.

This data is not included by default because then the SAPPHiRE package would
become to large. By running this script the data is added after installation.
This script downloads approximately 100 MB of data.

This can be run by simply running the installed script from the command line::

    $ extend_local_data

To make the script show information about what it will do add the help flag::

    $ extend_local_data --help

"""
import argparse

from .update_local_data import update_sublevel_tsv
from ..api import Network


def update_additional_local_tsv(progress=True):
    """Get location tsv data for all stations"""

    station_numbers = Network().station_numbers()

    for data_type in ['eventtime']:
        update_sublevel_tsv(data_type, station_numbers, progress)


def main():
    descr = """Add additional data to local data, or update already downloaded
             data. Making data available locally can greatly speed up a
             program which uses this data. Approximately 100 MB of data will
             be downloaded. The data contains the eventtime data, i.e. hourly
             number of events for all stations."""
    parser = argparse.ArgumentParser(description=descr)
    parser.parse_args()
    update_additional_local_tsv()
