import datetime
import os
import tempfile

from urllib.request import urlretrieve

import tables

from sapphire import esd

self_path = os.path.dirname(__file__)

test_data_path = os.path.join(self_path, 'test_data/esd_load_data.h5')
test_data_coincidences_path = os.path.join(self_path, 'test_data/esd_coincidence_data.h5')

events_source = os.path.join(self_path, 'test_data/events-s501-20120101.tsv')
weather_source = os.path.join(self_path, 'test_data/weather-s501-20120101.tsv')
singles_source = os.path.join(self_path, 'test_data/singles-s501-20170101.tsv')
lightning_source = os.path.join(self_path, 'test_data/lightning-knmi-20150717.tsv')
coincidences_source = os.path.join(self_path, 'test_data/coincidences-20160310.tsv')


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
        esd.load_data(datafile, '/', singles_source, 'singles')
        esd.load_data(datafile, '/', lightning_source, 'lightning')


def perform_esd_download_data(filename):
    """Load data from esd/api to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    start = datetime.datetime(2012, 1, 1, 0, 0, 0)
    end = datetime.datetime(2012, 1, 1, 0, 1, 0)
    singles_start = datetime.datetime(2017, 1, 1, 0, 0, 0)
    singles_end = datetime.datetime(2017, 1, 1, 0, 10, 0)
    lightning_start = datetime.datetime(2015, 7, 17, 0, 0, 0)
    lightning_end = datetime.datetime(2015, 7, 17, 0, 10, 0)

    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.download_data(datafile, '/', 501, start, end, type='events', progress=False)
        esd.download_data(datafile, '/', 501, start, end, type='weather', progress=False)
        esd.download_data(datafile, '/', 501, singles_start, singles_end, type='singles', progress=False)
        esd.download_lightning(datafile, '/', 4, lightning_start, lightning_end, progress=False)


def perform_load_coincidences(filename):
    """Load coincidence data from tsv (source) to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.load_coincidences(datafile, coincidences_source)


def perform_download_coincidences(filename):
    """Load data from esd/api to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    start = datetime.datetime(2016, 3, 10, 0, 0, 0)
    end = datetime.datetime(2016, 3, 10, 0, 1, 0)

    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.download_coincidences(datafile, stations=[501, 510], start=start, end=end, n=2, progress=False)


def create_and_store_test_data():
    """Create test data for future acceptance testing"""

    perform_esd_download_data(test_data_path)
    perform_download_coincidences(test_data_coincidences_path)
    urlretrieve(esd.get_weather_url().format(station_number=501, query='start=2012-01-01&end=2012-01-01+00:01:00'),
                weather_source)
    urlretrieve(esd.get_events_url().format(station_number=501, query='start=2012-01-01&end=2012-01-01+00:01:00'),
                events_source)
    urlretrieve(esd.get_singles_url().format(station_number=501, query='start=2017-01-01&end=2017-01-01+00:10:00'),
                singles_source)
    urlretrieve(esd.get_lightning_url().format(lightning_type=4, query='start=2015-07-17&end=2015-07-17+00:10:00'),
                lightning_source)
    urlretrieve(esd.get_coincidences_url().format(query='start=2016-03-10&end=2016-03-10+00:01:00&stations=501,+510&n=2'),
                coincidences_source)


if __name__ == '__main__':
    create_and_store_test_data()
