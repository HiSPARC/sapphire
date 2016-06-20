import os
import tempfile
from urllib import urlretrieve

import tables
import datetime

from sapphire import esd


self_path = os.path.dirname(__file__)
test_data_path = os.path.join(self_path, 'test_data/esd_load_data.h5')
test_data_coincidences_path = os.path.join(self_path,
                                           'test_data/esd_coincidence_data.h5')
events_source = os.path.join(self_path, 'test_data/events-s501-20120101.tsv')
weather_source = os.path.join(self_path, 'test_data/weather-s501-20120101.tsv')
coincidences_source = os.path.join(self_path,
                                   'test_data/coincidences-20160310.tsv')


def create_tempfile_path():
    """Create a temporary file, close it, and return the path"""

    f, path = tempfile.mkstemp(suffix='.h5')
    os.close(f)
    return path


def perform_load_data(filename):
    """Load data from tsv (source) to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.load_data(datafile, '/', events_source, 'events')
        esd.load_data(datafile, '/', weather_source, 'weather')


def perform_esd_download_data(filename):
    """Load data from esd/api to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    start = datetime.datetime(2012, 1, 1, 0, 0, 0)
    eind = datetime.datetime(2012, 1, 1, 0, 1, 0)

    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.download_data(datafile, '/', 501, start, eind, type='events',
                          progress=False)
        esd.download_data(datafile, '/', 501, start, eind, type='weather',
                          progress=False)


def perform_load_coincidences(filename):
    """Load coincidence data from tsv (source) to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.load_coincidences(datafile, coincidences_source)


def perform_download_coincidences(filename):
    """Load data from esd/api to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    start = datetime.datetime(2016, 3, 10, 0, 0, 0)
    eind = datetime.datetime(2016, 3, 10, 0, 1, 0)

    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.download_coincidences(datafile, stations=[501, 510],
                                  start=start, end=eind, n=2, progress=False)


def create_and_store_test_data():
    """Create test data for future acceptance testing"""

    perform_esd_download_data(test_data_path)
    perform_download_coincidences(test_data_coincidences_path)
    urlretrieve('http://data.hisparc.nl/data/501/weather/'
                '?download=True&start=2012-01-01&end=2012-01-01+00:01:00',
                weather_source)
    urlretrieve('http://data.hisparc.nl/data/501/events/'
                '?download=True&start=2012-01-01&end=2012-01-01+00:01:00',
                events_source)
    urlretrieve('http://data.hisparc.nl/data/network/coincidences/'
                '?download=True&start=2016-03-10&end=2016-03-10+00:01:00'
                '&stations=501,+510&n=2',
                coincidences_source)


if __name__ == '__main__':
    create_and_store_test_data()
