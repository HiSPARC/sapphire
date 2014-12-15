import unittest
import datetime

from mock import patch

from sapphire import time_util


class GPSTimeUtilTests(unittest.TestCase):
    def gpstime(self, gpstime):
        self.assertEqual(gpstime.gpstimestamp(), 1354320000)
        self.assertEqual(gpstime.description(), 'Sat Dec  1 00:00:00 2012')
        self.assertEqual(gpstime.datetime(), datetime.datetime(2012, 12, 1, 0, 0))

    def test_gpstime(self):
        self.gpstime(time_util.GPSTime(1354320000))
        self.gpstime(time_util.GPSTime(2012, 12, 1))
        self.gpstime(time_util.GPSTime(2012, 12, 1, 0, 0, 0))

    def test_incorrect_arguments(self):
        self.assertRaises(TypeError, time_util.GPSTime)
        self.assertRaises(TypeError, time_util.GPSTime, 2012, 12)

    @patch.object(time_util.GPSTime, 'description')
    def test_str_returns_description(self, mock_description):
        expected = "Foobar"
        mock_description.return_value = expected
        t = time_util.GPSTime(2014, 10, 27)
        actual = str(t)
        mock_description.assert_called_once_with()
        self.assertIs(actual, expected)


if __name__ == '__main__':
    unittest.main()
