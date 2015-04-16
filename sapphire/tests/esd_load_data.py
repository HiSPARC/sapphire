import os
import tempfile

import tables

from sapphire import esd


self_path = os.path.dirname(__file__)
test_data_path = os.path.join(self_path, 'test_data/esd_load_data.h5')
events_source = os.path.join(self_path, 'test_data/events-s501-20120101.csv')
weather_source = os.path.join(self_path, 'test_data/weather-s501-20120101.csv')


def create_tempfile_path():
    """Create a temporary file, close it, and return the path"""

    f, path = tempfile.mkstemp(suffix='.h5')
    os.close(f)
    return path


def perform_load_data(filename):
    """Load data from csv (source) to h5 (filename)"""

    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as datafile:
        esd.load_data(datafile, '/', events_source, 'events')
        esd.load_data(datafile, '/', weather_source, 'weather')


def create_and_store_test_data():
    """Create test data for future acceptance testing"""

    perform_load_data(test_data_path)


if __name__ == '__main__':
    create_and_store_test_data()
