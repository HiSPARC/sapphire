import os
import unittest

from mock import patch, ANY, sentinel, MagicMock
import tables

from sapphire import esd, api
from sapphire.tests.validate_results import validate_results

from esd_load_data import (create_tempfile_path,
                           test_data_path, test_data_coincidences_path,
                           perform_load_data, perform_load_coincidences,
                           perform_esd_download_data,
                           perform_download_coincidences)


class ESDTest(unittest.TestCase):

    def test_create_table(self):
        description = {'event_id': tables.UInt32Col(pos=0),
                       'timestamp': tables.Time32Col(pos=1),
                       'nanoseconds': tables.UInt32Col(pos=2),
                       'ext_timestamp': tables.UInt64Col(pos=3),
                       'pulseheights': tables.Int16Col(pos=4, shape=4),
                       'integrals': tables.Int32Col(pos=5, shape=4),
                       'n1': tables.Float32Col(pos=6),
                       'n2': tables.Float32Col(pos=7),
                       'n3': tables.Float32Col(pos=8),
                       'n4': tables.Float32Col(pos=9),
                       't1': tables.Float32Col(pos=10),
                       't2': tables.Float32Col(pos=11),
                       't3': tables.Float32Col(pos=12),
                       't4': tables.Float32Col(pos=13),
                       't_trigger': tables.Float32Col(pos=14)}
        file = MagicMock()
        result = esd._create_events_table(file, sentinel.group)
        file.create_table.assert_called_once_with(sentinel.group, 'events',
                                                  description,
                                                  createparents=True)
        self.assertEqual(result, file.create_table.return_value)

    def test_create_weather_table(self):
        description = {'event_id': tables.UInt32Col(pos=0),
                       'timestamp': tables.Time32Col(pos=1),
                       'temp_inside': tables.Float32Col(pos=2),
                       'temp_outside': tables.Float32Col(pos=3),
                       'humidity_inside': tables.Int16Col(pos=4),
                       'humidity_outside': tables.Int16Col(pos=5),
                       'barometer': tables.Float32Col(pos=6),
                       'wind_dir': tables.Int16Col(pos=7),
                       'wind_speed': tables.Int16Col(pos=8),
                       'solar_rad': tables.Int16Col(pos=9),
                       'uv': tables.Int16Col(pos=10),
                       'evapotranspiration': tables.Float32Col(pos=11),
                       'rain_rate': tables.Float32Col(pos=12),
                       'heat_index': tables.Int16Col(pos=13),
                       'dew_point': tables.Float32Col(pos=14),
                       'wind_chill': tables.Float32Col(pos=15)}
        file = MagicMock()
        result = esd._create_weather_table(file, sentinel.group)
        file.create_table.assert_called_once_with(sentinel.group, 'weather',
                                                  description,
                                                  createparents=True)
        self.assertEqual(result, file.create_table.return_value)

    def test__first_available_numbered_path(self):
        """Check if correct path is given if there is no existing h5."""

        self.assertEqual(esd._first_available_numbered_path(), 'data1.h5')
        # make data1.h5 and check if it returns data2.h5 then clean up..
        f = open('data1.h5', 'a')
        f.flush()
        f.close()
        self.assertEqual(esd._first_available_numbered_path(), 'data2.h5')
        os.remove('data1.h5')

    def test_unsupported_type(self):
        """Check for Exception for unsupported data types"""

        self.assertRaises(ValueError, esd.load_data, None, None, None, 'bad')
        self.assertRaises(ValueError, esd.download_data, None, None, 501,
                          type='bad')

    def test_start_end_values(self):
        """Check for RuntimeError for impossible end=value with start=None"""

        self.assertRaises(RuntimeError, esd.download_data, None, None, 501,
                          start=None, end='a_value')
        self.assertRaises(RuntimeError, esd.download_coincidences, None,
                          start=None, end="a_value")

    def test_load_data_output(self):
        """Load data tsv into hdf5 and verify the output"""

        output_path = create_tempfile_path()
        perform_load_data(output_path)
        validate_results(self, test_data_path, output_path)
        os.remove(output_path)

    def test_load_coincidences_output(self):
        """Load coincidences tsv into hdf5 and verify the output"""

        output_path = create_tempfile_path()
        perform_load_coincidences(output_path)
        validate_results(self, test_data_coincidences_path, output_path)
        os.remove(output_path)

    @patch.object(esd, 'download_data')
    @patch.object(tables, 'open_file')
    def test_quick_download(self, mock_open_file, mock_download_data):
        """Test esd.quick_download()"""

        esd.quick_download(501)
        mock_open_file.assert_called_once_with('data1.h5', 'w')
        mock_download_data.assert_called_once_with(ANY, None, 501, None)

    @unittest.skipUnless(api.API.check_connection(),
                         "Internet connection required")
    def test_download_data(self):
        """Download data and validate results"""

        output_path = create_tempfile_path()
        perform_esd_download_data(output_path)
        validate_results(self, test_data_path, output_path)
        os.remove(output_path)

    @unittest.skipUnless(api.API.check_connection(),
                         "Internet connection required")
    def test_download_coincidences(self):
        """Download coincidence data from esd and validate results"""

        output_path = create_tempfile_path()
        perform_download_coincidences(output_path)
        validate_results(self, test_data_coincidences_path, output_path)
        os.remove(output_path)


if __name__ == '__main__':
    unittest.main()
