import unittest

from mock import patch, sentinel
from numpy import isnan, nan, array, all
from datetime import date, datetime

from sapphire.analysis import calibration
from sapphire.transformations.clock import datetime_to_gps


class DetectorTimingTests(unittest.TestCase):

    @patch.object(calibration, 'fit_timing_offset')
    def test_determine_detector_timing_offset(self, mock_fit):
        # Empty list
        offset = calibration.determine_detector_timing_offset([])
        self.assertTrue(all(isnan(offset)))

        dt = array([-10, 0, 10])

        # Good result
        mock_fit.return_value = (1., 2.)
        offset, _ = calibration.determine_detector_timing_offset(dt)
        self.assertEqual(offset, 1.)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=.6)
        self.assertEqual(offset, 3.)

        mock_fit.return_value = (-1.5, 5.)
        offset, _ = calibration.determine_detector_timing_offset(dt)
        self.assertEqual(offset, -1.5)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=.6)

        self.assertEqual(offset, 0.5)

        mock_fit.return_value = (250., 100.)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=.6)
        self.assertTrue(isnan(offset))
        mock_fit.return_value = (-150., 100.)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=.6)
        self.assertTrue(isnan(offset))

        mock_fit.return_value = (nan, nan)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=.6)
        self.assertTrue(isnan(offset))


class StationTimingTests(unittest.TestCase):

    @patch.object(calibration, 'percentile')
    @patch.object(calibration, 'fit_timing_offset')
    def test_determine_station_timing_offset(self, mock_fit, mock_percentile):
        mock_percentile.return_value = (-50., 50.)

        # Empty list
        offset = calibration.determine_station_timing_offset([])
        self.assertTrue(all(isnan(offset)))

        # Good result
        mock_fit.return_value = (1., 5.)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertEqual(offset, 1.)
        mock_percentile.assert_called_once_with([sentinel.dt], [0.5, 99.5])
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt],
                                                                dz=.6)
        self.assertEqual(offset, 3.)

        mock_fit.return_value = (-1.5, 5.)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertEqual(offset, -1.5)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt],
                                                                dz=.6)
        self.assertEqual(offset, 0.5)

        mock_fit.return_value = (2500., 100.)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))
        mock_fit.return_value = (-1500., 100.)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))

        mock_fit.return_value = (nan, nan)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))


class BestReferenceTests(unittest.TestCase):

    def test_determine_best_reference(self):
        # Tie
        filters = array([[True, True, False], [True, False, True],
                         [False, True, True], [True, True, False]])
        self.assertEqual(calibration.determine_best_reference(filters), 0)

        # 1 has most matches
        filters = array([[True, False, False], [True, True, True],
                         [False, False, False], [True, True, False]])
        self.assertEqual(calibration.determine_best_reference(filters), 1)

        # Another winner
        filters = array([[True, True, False], [True, False, True],
                         [False, True, True], [True, True, True]])
        self.assertEqual(calibration.determine_best_reference(filters), 3)

        # Not yet support number of detectors
        filters = array([[True, True, False], [True, False, True]])
        self.assertRaises(IndexError, calibration.determine_best_reference,
                          filters)


class SplitDatetimeRangeTests(unittest.TestCase):

    def test_split_range(self):
        # 101 days
        start = date(2016, 1, 1)
        end_5days = date(2016, 1, 6)
        end_100days = date(2016, 4, 11)

        # no step, dates:
        result = list(calibration.datetime_range(start, end_5days))
        self.assertEqual(len(result), 5)
        begin, _ = result[0]
        _, end = result[-1]
        self.assertEqual(begin, start)
        self.assertEqual(end, end_5days)

        # single interval
        result = list(calibration.datetime_range(start, end_5days, 5))
        self.assertEqual(len(result), 1)
        begin, end = result[0]
        self.assertEqual(begin, start)
        self.assertEqual(end, end_5days)

        # split an even interval in two parts
        result = list(calibration.datetime_range(start, end_5days, 2))
        self.assertEqual(len(result), 2)
        begin, _ = result[0]
        _, end = result[-1]
        self.assertEqual(begin, start)
        self.assertEqual(end, end_5days)

        # split large interval, remainder = 0
        result = list(calibration.datetime_range(start, end_100days, 10))
        self.assertEqual(len(result), 10)
        begin, _ = result[0]
        _, end = result[-1]
        self.assertEqual(begin, start)
        self.assertEqual(end, end_100days)

        # split large interval, divide remainder
        result = list(calibration.datetime_range(start, end_100days, 7))
        self.assertEqual(len(result), 14)
        begin, _ = result[0]
        _, end = result[-1]
        self.assertEqual(begin, start)
        self.assertEqual(end, end_100days)


class DetermineTimingOffsetsTests(unittest.TestCase):

    def test_intervals(self):
        first = ('%d\t51.0\t4.0\t0.0\n' % datetime_to_gps(datetime(2014, 1, 1)) +
                 '%d\t51.0\t4.0\t0.0\n' % datetime_to_gps(datetime(2014, 1, 21)) +
                 '%d\t51.0\t4.0\t0.0\n' % datetime_to_gps(datetime(2014, 3, 1)))

        second = ('%d\t51.0\t4.0\t0.0\n' % datetime_to_gps(datetime(2014, 1, 7)) +
                  '%d\t51.0\t4.0\t0.0\n' % datetime_to_gps(datetime(2014, 1, 21)) +
                  '%d\t51.0\t4.0\t0.0\n' % datetime_to_gps(datetime(2014, 3, 5)))
        # WIP
        first, second = second, first


if __name__ == '__main__':
    unittest.main()
