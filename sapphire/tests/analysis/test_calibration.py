import unittest

from mock import patch, sentinel
from numpy import isnan, nan, array

from sapphire.analysis import calibration


class DetectorTimingTests(unittest.TestCase):

    @patch.object(calibration, 'fit_timing_offset')
    def test_determine_detector_timing_offset(self, mock_fit):
        # Empty list
        offset = calibration.determine_detector_timing_offset([])
        self.assertTrue(isnan(offset))

        # Good result
        mock_fit.return_value = 1.
        offset = calibration.determine_detector_timing_offset([sentinel.dt])
        self.assertEqual(offset, 1.)
        offset = calibration.determine_detector_timing_offset([sentinel.dt],
                                                              dz=.6)
        self.assertEqual(offset, 3.)

        mock_fit.return_value = -1.5
        offset = calibration.determine_detector_timing_offset([sentinel.dt])
        self.assertEqual(offset, -1.5)
        offset = calibration.determine_detector_timing_offset([sentinel.dt],
                                                              dz=.6)
        self.assertEqual(offset, 0.5)

        mock_fit.return_value = 250.
        offset = calibration.determine_detector_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))
        mock_fit.return_value = -150.
        offset = calibration.determine_detector_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))

        mock_fit.return_value = nan
        offset = calibration.determine_detector_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))


class StationTimingTests(unittest.TestCase):

    @patch.object(calibration, 'percentile')
    @patch.object(calibration, 'fit_timing_offset')
    def test_determine_station_timing_offset(self, mock_fit, mock_percentile):
        mock_percentile.return_value = (-50., 50.)

        # Empty list
        offset = calibration.determine_station_timing_offset([])
        self.assertTrue(isnan(offset))

        # Good result
        mock_fit.return_value = 1.
        offset = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertEqual(offset, 1.)
        mock_percentile.assert_called_once_with([sentinel.dt], [2, 98])
        offset = calibration.determine_station_timing_offset([sentinel.dt],
                                                             dz=.6)
        self.assertEqual(offset, 3.)

        mock_fit.return_value = -1.5
        offset = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertEqual(offset, -1.5)
        offset = calibration.determine_station_timing_offset([sentinel.dt],
                                                             dz=.6)
        self.assertEqual(offset, 0.5)

        mock_fit.return_value = 2500.
        offset = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))
        mock_fit.return_value = -1500.
        offset = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertTrue(isnan(offset))

        mock_fit.return_value = nan
        offset = calibration.determine_station_timing_offset([sentinel.dt])
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


if __name__ == '__main__':
    unittest.main()
