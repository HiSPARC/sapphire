import os
import unittest

from mock import sentinel, MagicMock
import tables

from sapphire import esd

from esd_load_data import (create_tempfile_path, perform_load_data,
                           test_data_path)


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

    def test_esd_output(self):
        """Perform a simulation and verify the output"""

        output_path = create_tempfile_path()
        perform_load_data(output_path)
        self.validate_results(test_data_path, output_path)
        os.remove(output_path)

    def validate_results(self, expected_path, actual_path):
        """Validate simulation results"""

        with tables.open_file(expected_path) as expected_file:
            with tables.open_file(actual_path) as actual_file:
                for table in ('/events', '/weather'):
                    self.validate_table(table, expected_file, actual_file)

    def validate_table(self, table, expected_file, actual_file):
        """Verify that two tables are identical"""

        expected_node = expected_file.get_node(table)
        actual_node = actual_file.get_node(table)

        for colname in expected_node.colnames:
            expected_col = expected_node.col(colname)
            actual_col = actual_node.col(colname)
            if expected_col.shape == actual_col.shape:
                self.assertTrue((expected_col == actual_col).all())
            else:
                self.fail("Columns do not have the same length.")


if __name__ == '__main__':
    unittest.main()
