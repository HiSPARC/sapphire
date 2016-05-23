import unittest

from mock import patch, sentinel, MagicMock, Mock, call
from numpy import isnan, nan, array, all
from datetime import datetime, date

import sapphire
from sapphire.analysis import calibration
from sapphire.transformations.clock import datetime_to_gps
from sapphire.utils import c


class DetectorTimingTests(unittest.TestCase):

    @patch.object(calibration, 'fit_timing_offset')
    def test_determine_detector_timing_offset(self, mock_fit):
        # Empty list
        offset = calibration.determine_detector_timing_offset([])
        self.assertTrue(all(isnan(offset)))

        dt = array([-10, 0, 10])
        dz = 0.6
        dzc = dz / c

        # Good result
        mock_fit.return_value = (1., 2.)
        offset, _ = calibration.determine_detector_timing_offset(dt)
        self.assertEqual(offset, 1.)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=dz)
        self.assertEqual(offset, 1. + dzc)

        mock_fit.return_value = (-1.5, 5.)
        offset, _ = calibration.determine_detector_timing_offset(dt)
        self.assertEqual(offset, -1.5)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=dz)
        self.assertEqual(offset, -1.5 + dzc)

        mock_fit.return_value = (250., 100.)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=dz)
        self.assertTrue(isnan(offset))
        mock_fit.return_value = (-150., 100.)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=dz)
        self.assertTrue(isnan(offset))

        mock_fit.return_value = (nan, nan)
        offset, _ = calibration.determine_detector_timing_offset(dt, dz=dz)
        self.assertTrue(isnan(offset))


class StationTimingTests(unittest.TestCase):

    @patch.object(calibration, 'percentile')
    @patch.object(calibration, 'fit_timing_offset')
    def test_determine_station_timing_offset(self, mock_fit, mock_percentile):
        mock_percentile.return_value = (-50., 50.)
        dz = 0.6
        dzc = dz / c

        # Empty list
        offset = calibration.determine_station_timing_offset([])
        self.assertTrue(all(isnan(offset)))

        # Good result
        mock_fit.return_value = (1., 5.)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertEqual(offset, 1.)
        mock_percentile.assert_called_once_with([sentinel.dt], [0.5, 99.5])
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt],
                                                                dz=dz)
        self.assertEqual(offset, 1. + dzc)

        mock_fit.return_value = (-1.5, 5.)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt])
        self.assertEqual(offset, -1.5)
        offset, _ = calibration.determine_station_timing_offset([sentinel.dt],
                                                                dz=dz)
        self.assertEqual(offset, -1.5 + dzc)

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

    def test_pairwise(self):
        result = list(calibration.pairwise([1, 2, 3, 4]))
        self.assertEqual(result, [(1, 2), (2, 3), (3, 4)])


class DetermineStationTimingOffsetsTests(unittest.TestCase):

    def setUp(self):
        stations = [501, 102, 105, 8001]
        self.off = calibration.DetermineStationTimingOffsets(stations=stations, data=sentinel.data,
                                                             progress=sentinel.progress, force_stale=True)

    def test_init(self):
        self.assertEqual(self.off.progress, sentinel.progress)
        self.assertEqual(self.off.data, sentinel.data)
        self.assertIsInstance(self.off.cluster, sapphire.clusters.HiSPARCStations)

    def test_init_network(self):
        off = calibration.DetermineStationTimingOffsets(force_stale=True)
        self.assertIsInstance(off.cluster, sapphire.clusters.HiSPARCNetwork)

    def test_read_dt(self):
        self.off.data = MagicMock()
        table_mock = MagicMock()
        self.off.data.get_node.return_value = table_mock
        station = 502
        ref_station = 501
        start = datetime(2014, 1, 1)
        end = datetime(2016, 12, 31)
        self.off.read_dt(station, ref_station, start, end)
        table_path = '/time_deltas/station_%d/station_%d' % (ref_station,
                                                             station)
        table_name = 'time_deltas'
        self.off.data.get_node.assert_called_once_with(table_path, table_name)
        self.assertTrue(table_mock.read_where.called)

    def test_station_pairs_within_max_distance(self):
        results = list(self.off.get_station_pairs_within_max_distance())
        self.assertEqual([(102, 105)], results)
        # force unsorted order of stations:

    def test_station_pairs_wrong_order(self):
        stations = [105, 102, 8001, 501]
        self.off = calibration.DetermineStationTimingOffsets(stations=stations, data=sentinel.data,
                                                             progress=sentinel.progress, force_stale=True)
        results = list(self.off.get_station_pairs_within_max_distance())
        self.assertEqual([(102, 105)], results)

    @patch.object(calibration, 'Station')
    def test_get_gps_timestamps(self, mock_station):
        self.off._get_gps_timestamps(sentinel.station)
        mock_station.assert_called_once_with(sentinel.station, force_stale=self.off.force_stale)

    @patch.object(calibration, 'Station')
    def test_get_electronics_timestamp(self, mock_station):
        self.off._get_electronics_timestamps(sentinel.station)
        mock_station.assert_called_once_with(sentinel.station, force_stale=self.off.force_stale)

    def test_get_r_dz(self):
        r_102_105 = 82.2162521322  # 2014,1,1
        dz_102_105 = -4.13568408095
        r, dz = self.off._get_r_dz(datetime(2014, 1, 1).date(), 102, 105)
        self.assertAlmostEqual(r, r_102_105)
        self.assertAlmostEqual(dz, dz_102_105)

    def test_determine_interval(self):
        combinations = ((0., 7),
                        (50., 7),
                        (200., 37),
                        (1000., 229))
        for r, ref_int in combinations:
            self.assertEqual(self.off._determine_interval(r), ref_int)

    def test_get_cuts(self):
        gps_station = (datetime_to_gps(datetime(2014, 1, 1, 10, 3)),
                       datetime_to_gps(datetime(2014, 3, 1, 11, 32)))
        gps_ref_station = (datetime_to_gps(datetime(2014, 1, 5, 0, 1, 1)),
                           datetime_to_gps(datetime(2014, 3, 5, 3, 34, 4)))
        elec_station = (datetime_to_gps(datetime(2014, 1, 3, 3, 34, 3)),
                        datetime_to_gps(datetime(2014, 3, 5, 23, 59, 59)))
        elec_ref_station = (datetime_to_gps(datetime(2014, 1, 9, 0, 0, 0)),
                            datetime_to_gps(datetime(2014, 3, 15, 1, 2, 3)))
        gps_mock = Mock()
        elec_mock = Mock()

        gps_mock.side_effect = [array(gps_station), array(gps_ref_station)]
        elec_mock.side_effect = [array(elec_station), array(elec_ref_station)]

        self.off._get_electronics_timestamps = elec_mock
        self.off._get_gps_timestamps = gps_mock

        cuts = self.off._get_cuts(sentinel.station, sentinel.ref_station)

        elec_mock.assert_has_calls([call(sentinel.ref_station), call(sentinel.station)], any_order=True)
        gps_mock.assert_has_calls([call(sentinel.ref_station), call(sentinel.station)], any_order=True)

        self.assertEqual(len(cuts), 9)
        self.assertItemsEqual(sorted(cuts), cuts)
        self.assertEqual(cuts[0], datetime(2014, 1, 1))
        today = datetime.now()
        self.assertEqual(cuts[-1], datetime(today.year, today.month, today.day))

    def test_get_left_and_right_bounds(self):
        cuts = (datetime(2014, 1, 1),
                datetime(2015, 1, 1),
                datetime(2015, 1, 5),
                datetime(2015, 1, 10))
        combinations = [(datetime(2015, 1, 1), 7, datetime(2015, 1, 1), datetime(2015, 1, 4)),
                        (datetime(2015, 1, 3), 7, datetime(2015, 1, 1), datetime(2015, 1, 4)),
                        (datetime(2015, 1, 3).date(), 7, datetime(2015, 1, 1), datetime(2015, 1, 4)),
                        (datetime(2015, 1, 5), 7, datetime(2015, 1, 5), datetime(2015, 1, 9)),
                        (datetime(2015, 1, 10), 7, datetime(2015, 1, 5), datetime(2015, 1, 10))]
        for d, days, ref_left, ref_right in combinations:
            left, right = self.off._get_left_and_right_bounds(cuts, d, days)
            self.assertEqual(left, ref_left)
            self.assertEqual(right, ref_right)

    def test_determine_first_and_last_date(self):
        date = datetime(2015, 1, 2)
        self.off._determine_interval = Mock()
        self.off._get_r_dz = Mock()
        self.off._get_cuts = Mock()
        self.off._get_left_and_right_bounds = Mock()

        self.off._get_r_dz.return_value = sentinel.r, sentinel.dz
        self.off._get_cuts.return_value = sentinel.cuts
        self.off._determine_interval.return_value = sentinel.interval

        self.off.determine_first_and_last_date(date, sentinel.station, sentinel.ref_station)

        self.off._get_cuts.assert_called_once_with(sentinel.station, sentinel.ref_station)
        self.off._get_r_dz.assert_called_once_with(date, sentinel.station, sentinel.ref_station)
        self.off._determine_interval.assert_called_once_with(sentinel.r)
        self.off._get_left_and_right_bounds.assert_called_once_with(sentinel.cuts, date, sentinel.interval)

        # Test if datetime.date objects are handled correctly
        date = datetime(2015, 1, 2).date()
        self.off.determine_first_and_last_date(date, sentinel.station, sentinel.ref_station)

    def test_datetime(self):
        d = datetime(2000, 1, 2)
        self.assertEqual(datetime(2000, 1, 2), self.off._datetime(d))

        d = datetime(2000, 1, 2).date()
        self.assertEqual(datetime(2000, 1, 2), self.off._datetime(d))

        d = datetime(2000, 1, 2, 3, 4, 5, 6)
        self.assertEqual(datetime(2000, 1, 2), self.off._datetime(d))

        d = datetime(2000, 1, 2, 3, 4, 5, 6).date()
        self.assertEqual(datetime(2000, 1, 2), self.off._datetime(d))

    @patch.object(calibration, 'determine_station_timing_offset')
    def test_determine_station_timing_offset(self, mock_det_offset):
        date = datetime(2015, 1, 2)
        self.off._get_r_dz = Mock()
        self.off.determine_first_and_last_date = Mock()
        self.off.read_dt = Mock()

        self.off._get_r_dz.return_value = sentinel.r, sentinel.dz
        self.off.determine_first_and_last_date.return_value = (0, 0)
        self.off.read_dt.return_value = 1000 * [0.]
        mock_det_offset.return_value = (10., 1.)

        offsets = self.off.determine_station_timing_offset(date, sentinel.station, sentinel.ref_station)

        self.off._get_r_dz.assert_called_once_with(date, sentinel.station, sentinel.ref_station)
        self.off.determine_first_and_last_date.assert_called_once_with(date, sentinel.station, sentinel.ref_station)

        self.assertEqual(offsets, (10., 1.))

        self.off.read_dt.return_value = 90 * [0.]
        offsets = self.off.determine_station_timing_offset(date, sentinel.station, sentinel.ref_station)
        self.assertEqual(offsets, (nan, nan))


if __name__ == '__main__':
    unittest.main()
