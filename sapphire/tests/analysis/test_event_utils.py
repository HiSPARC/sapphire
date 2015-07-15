import unittest
import warnings

from mock import MagicMock, patch, sentinel
from numpy import isnan, nan

from sapphire.analysis import event_utils


class StationDensityTests(unittest.TestCase):

    @patch.object(event_utils, 'detector_densities')
    @patch.object(event_utils, 'get_detector_ids')
    def test_station_density(self, mock_detector_ids, mock_detector_densities):
        mock_detector_ids.return_value = range(4)
        mock_detector_densities.return_value = [1, 1, 2, 2]
        self.assertEqual(event_utils.station_density(sentinel.event, range(4), sentinel.station), 1.5)
        mock_detector_densities.return_value = [1, 1, nan, nan]
        self.assertEqual(event_utils.station_density(sentinel.event, range(4), sentinel.station), 1)
        self.assertEqual(mock_detector_ids.call_count, 0)
        self.assertEqual(event_utils.station_density(sentinel.event), 1)
        mock_detector_ids.assert_called_once_with(None, sentinel.event)
        self.assertEqual(event_utils.station_density(sentinel.event, station=sentinel.station), 1)
        mock_detector_ids.assert_called_with(sentinel.station, sentinel.event)


class DetectorDensitiesTests(unittest.TestCase):

    @patch.object(event_utils, 'detector_density')
    @patch.object(event_utils, 'get_detector_ids')
    def test_detector_densities(self, mock_detector_ids, mock_detector_density):
        mock_detector_ids.return_value = range(4)
        mock_detector_density.return_value = sentinel.density
        self.assertEqual(event_utils.detector_densities(sentinel.event, range(4)),
                         [sentinel.density] * 4)
        self.assertEqual(event_utils.detector_densities(sentinel.event, range(2)),
                         [sentinel.density] * 2)
        self.assertEqual(mock_detector_ids.call_count, 0)
        self.assertEqual(event_utils.detector_densities(sentinel.event),
                         [sentinel.density] * 4)
        mock_detector_ids.assert_called_once_with(None, sentinel.event)
        self.assertEqual(event_utils.detector_densities(sentinel.event, station=sentinel.station),
                         [sentinel.density] * 4)
        mock_detector_ids.assert_called_with(sentinel.station, sentinel.event)


class DetectorDensityTests(unittest.TestCase):

    def setUp(self):
        self.event = MagicMock()
        self.station = MagicMock()

    def test_detector_density(self):
        self.event.__getitem__.side_effect = lambda name: 2
        self.assertEqual(event_utils.detector_density(self.event, 0), 4)
        self.event.__getitem__.assert_called_with('n1')

    def test_no_good_detector_density(self):
        self.event.__getitem__.side_effect = lambda name: -999
        self.assertTrue(isnan(event_utils.detector_density(self.event, 0)))
        self.event.__getitem__.assert_called_with('n1')


class StationArrivalTimeTests(unittest.TestCase):

    @patch.object(event_utils, 'detector_arrival_times')
    @patch.object(event_utils, 'get_detector_ids')
    def test_station_arrival_time(self, mock_detector_ids, mock_detector_arrival_times):
        mock_detector_ids.return_value = range(4)
        mock_detector_arrival_times.return_value = [7.5, 5., 2.5, 5.]
        event_dict = {'t_trigger': 10, 'ext_timestamp': 1000}
        event = MagicMock()
        event.__getitem__.side_effect = lambda name: event_dict[name]
        ref_ets = 500
        rel_arrival_time = event_dict['ext_timestamp'] - ref_ets - event_dict['t_trigger']

        self.assertEqual(event_utils.station_arrival_time(event, ref_ets, range(4), sentinel.offsets, sentinel.station),
                         rel_arrival_time + 2.5)
        self.assertEqual(mock_detector_ids.call_count, 0)
        self.assertEqual(event_utils.station_arrival_time(event, ref_ets, None, sentinel.offsets),
                         rel_arrival_time + 2.5)
        mock_detector_ids.assert_called_once_with(None, event)
        self.assertEqual(event_utils.station_arrival_time(event, ref_ets, None, sentinel.offsets, sentinel.station),
                         rel_arrival_time + 2.5)
        mock_detector_ids.assert_called_with(sentinel.station, event)

    @patch.object(event_utils, 'detector_arrival_times')
    @patch.object(event_utils, 'get_detector_ids')
    def test_nan_station_arrival_time(self, mock_detector_ids, mock_detector_arrival_times):
        mock_detector_ids.return_value = range(4)
        mock_detector_arrival_times.return_value = [7.5, 5., nan, 5.]
        event_dict = {'t_trigger': 10, 'ext_timestamp': 1000}
        event = MagicMock()
        event.__getitem__.side_effect = lambda name: event_dict[name]
        ref_ets = 500
        rel_arrival_time = event_dict['ext_timestamp'] - ref_ets - event_dict['t_trigger']

        self.assertEqual(event_utils.station_arrival_time(event, ref_ets, None, sentinel.offsets, sentinel.station),
                         rel_arrival_time + 5)
        event_dict['t_trigger'] = -999
        self.assertTrue(isnan(event_utils.station_arrival_time(event, ref_ets, None, sentinel.offsets, sentinel.station)))
        event_dict['t_trigger'] = nan
        self.assertTrue(isnan(event_utils.station_arrival_time(event, ref_ets, None, sentinel.offsets, sentinel.station)))
        event_dict['t_trigger'] = 10
        mock_detector_arrival_times.return_value = [nan, nan, nan, nan]
        with warnings.catch_warnings(record=True) as warned:
            self.assertTrue(isnan(event_utils.station_arrival_time(event, ref_ets, None, sentinel.offsets, sentinel.station)))
        self.assertEqual(len(warned), 1)


class RelativeDetectorArrivalTimesTests(unittest.TestCase):

    @patch.object(event_utils, 'detector_arrival_times')
    @patch.object(event_utils, 'get_detector_ids')
    def test_relative_detector_arrival_times(self, mock_detector_ids, mock_detector_arrival_times):
        mock_detector_ids.return_value = range(4)
        mock_detector_arrival_times.return_value = [7.5, 5., 2.5, 5.]
        event_dict = {'t_trigger': 10, 'ext_timestamp': 1000}
        event = MagicMock()
        event.__getitem__.side_effect = lambda name: event_dict[name]
        ref_ets = 500
        rel_arrival_time = event_dict['ext_timestamp'] - ref_ets - event_dict['t_trigger']

        self.assertEqual(event_utils.relative_detector_arrival_times(event, 500, range(4), sentinel.offsets, sentinel.station),
                         [rel_arrival_time + t for t in mock_detector_arrival_times()])
        self.assertEqual(mock_detector_ids.call_count, 0)
        self.assertEqual(event_utils.relative_detector_arrival_times(event, 500, None, sentinel.offsets),
                         [rel_arrival_time + t for t in mock_detector_arrival_times()])
        mock_detector_ids.assert_called_once_with(None, event)
        self.assertEqual(event_utils.relative_detector_arrival_times(event, 500, None, sentinel.offsets, sentinel.station),
                         [rel_arrival_time + t for t in mock_detector_arrival_times()])
        mock_detector_ids.assert_called_with(sentinel.station, event)

    @patch.object(event_utils, 'detector_arrival_times')
    @patch.object(event_utils, 'get_detector_ids')
    def test_nan_relative_detector_arrival_times(self, mock_detector_ids, mock_detector_arrival_times):
        mock_detector_ids.return_value = range(4)
        mock_detector_arrival_times.return_value = [7.5, 5., 5., nan]
        event_dict = {'t_trigger': 10, 'ext_timestamp': 1000}
        event = MagicMock()
        event.__getitem__.side_effect = lambda name: event_dict[name]
        ref_ets = 500
        rel_arrival_time = event_dict['ext_timestamp'] - ref_ets - event_dict['t_trigger']

        result = event_utils.relative_detector_arrival_times(event, 500, None, sentinel.offsets, sentinel.station)
        self.assertEqual(result[:-1], [rel_arrival_time + t for t in mock_detector_arrival_times()[:-1]])
        self.assertTrue(isnan(result[-1]))
        event_dict['t_trigger'] = -999
        self.assertTrue(isnan(event_utils.relative_detector_arrival_times(event, 500, None, sentinel.offsets, sentinel.station)).all())
        event_dict['t_trigger'] = nan
        self.assertTrue(isnan(event_utils.relative_detector_arrival_times(event, 500, None, sentinel.offsets, sentinel.station)).all())
        event_dict['t_trigger'] = 10
        mock_detector_arrival_times.return_value = [nan, nan, nan, nan]
        self.assertTrue(isnan(event_utils.relative_detector_arrival_times(event, 500, None, sentinel.offsets, sentinel.station)).all())


class DetectorArrivalTimesTests(unittest.TestCase):

    @patch.object(event_utils, 'detector_arrival_time')
    @patch.object(event_utils, 'get_detector_ids')
    def test_detector_arrival_times(self, mock_detector_ids, mock_detector_arrival_time):
        mock_detector_ids.return_value = range(4)

        mock_detector_arrival_time.return_value = sentinel.time
        self.assertEqual(event_utils.detector_arrival_times(sentinel.event, range(4)),
                         [sentinel.time] * 4)
        self.assertEqual(event_utils.detector_arrival_times(sentinel.event, range(2)),
                         [sentinel.time] * 2)
        self.assertEqual(mock_detector_ids.call_count, 0)
        self.assertEqual(event_utils.detector_arrival_times(sentinel.event),
                         [sentinel.time] * 4)
        mock_detector_ids.assert_called_once_with(None, sentinel.event)
        self.assertEqual(event_utils.detector_arrival_times(sentinel.event, station=sentinel.station),
                         [sentinel.time] * 4)
        mock_detector_ids.assert_called_with(sentinel.station, sentinel.event)


class DetectorArrivalTimeTests(unittest.TestCase):

    def setUp(self):
        self.event = MagicMock()
        self.offsets = [1, 2, 3, 4]

    def test_detector_arrival_time(self):
        self.event.__getitem__.side_effect = lambda name: 2.5
        self.assertEqual(event_utils.detector_arrival_time(self.event, 0), 2.5)
        self.event.__getitem__.assert_called_with('t1')
        self.assertEqual(event_utils.detector_arrival_time(self.event, 0, self.offsets), 1.5)
        self.event.__getitem__.assert_called_with('t1')
        self.assertEqual(event_utils.detector_arrival_time(self.event, 1, self.offsets), 0.5)
        self.event.__getitem__.assert_called_with('t2')

    def test_no_good_detector_arrival_time(self):
        self.event.__getitem__.side_effect = lambda name: -999
        self.assertTrue(isnan(event_utils.detector_arrival_time(self.event, 0)))
        self.event.__getitem__.assert_called_with('t1')
        self.assertTrue(isnan(event_utils.detector_arrival_time(self.event, 0, self.offsets)))
        self.event.__getitem__.assert_called_with('t1')


class GetDetectorIdsTests(unittest.TestCase):

    def test_get_detector_ids(self):
        self.assertEqual(event_utils.get_detector_ids(), range(4))
        station = MagicMock()
        station.detectors.__len__.return_value = 2
        self.assertEqual(event_utils.get_detector_ids(station=station), range(2))
        event = MagicMock()
        event.__getitem__.side_effect = lambda name: [10, 100, 40, -1]
        self.assertEqual(event_utils.get_detector_ids(event=event), range(3))
        self.assertEqual(event_utils.get_detector_ids(station=station, event=event), range(2))


if __name__ == '__main__':
    unittest.main()
