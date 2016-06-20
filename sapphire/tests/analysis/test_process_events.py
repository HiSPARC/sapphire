import unittest
import tempfile
import os
import shutil
import warnings
import operator

import tables
from mock import Mock

from sapphire.analysis import process_events


TEST_DATA_FILE = 'test_data/process_events.h5'
DATA_GROUP = '/s501'


class ProcessEventsTests(unittest.TestCase):
    def setUp(self):
        warnings.filterwarnings('ignore')
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.open_file(self.data_path, 'a')
        self.proc = process_events.ProcessEvents(self.data, DATA_GROUP, progress=False)

    def tearDown(self):
        warnings.resetwarnings()
        self.data.close()
        os.remove(self.data_path)

    def test_get_traces_for_event(self):
        event = self.proc.source[0]
        self.assertEqual(self.proc.get_traces_for_event(event)[12][3], 1334)

    def test__find_unique_row_ids(self):
        ext_timestamps = self.proc.source.col('ext_timestamp')
        enumerated_timestamps = list(enumerate(ext_timestamps))
        enumerated_timestamps.sort(key=operator.itemgetter(1))
        ids_in = [id for id, _ in enumerated_timestamps]
        ids = self.proc._find_unique_row_ids(enumerated_timestamps)
        self.assertEqual(ids, ids_in)

        enumerated_timestamps = [(0, 1), (1, 1), (3, 2), (2, 2)]
        ids = self.proc._find_unique_row_ids(enumerated_timestamps)
        self.assertEqual(ids, [0, 3])

        # Must be sorted by timestamp or the result will be differenct.
        enumerated_timestamps = [(0, 1), (3, 2), (1, 1), (2, 2)]
        ids = self.proc._find_unique_row_ids(enumerated_timestamps)
        self.assertNotEqual(ids, [0, 3])

    def test__reconstruct_time_from_traces(self):
        event = self.proc.source[10]
        times = self.proc._reconstruct_time_from_traces(event)
        self.assertEqual(times[0], 162.5)
        self.assertEqual(times[2], -999)

        event['pulseheights'][0] = -1
        times = self.proc._reconstruct_time_from_traces(event)
        self.assertEqual(times[0], -1)

    def test__reconstruct_time_from_trace(self):
        trace = [220, 222, 224, 222, 220]
        self.assertEqual(self.proc._reconstruct_time_from_trace(trace, 200), 0)
        self.assertEqual(self.proc._reconstruct_time_from_trace(trace, 203), 2)
        self.assertEqual(self.proc._reconstruct_time_from_trace(trace, 205), -999)

    def test_first_above_threshold(self):
        trace = [0, 2, 4, 2, 0]
        self.assertEqual(self.proc.first_above_threshold(trace, 1), 1)
        self.assertEqual(self.proc.first_above_threshold(trace, 3), 2)
        self.assertEqual(self.proc.first_above_threshold(trace, 4), 2)
        self.assertEqual(self.proc.first_above_threshold(trace, 5), -999)

#     @patch.object(process_events.FindMostProbableValueInSpectrum, 'find_mpv')
    def test__process_pulseintegrals(self):
        self.proc.limit = 1
#         mock_find_mpv.return_value = (-999, False)
        # Because of small data sample fit fails for detector 1
        self.assertEqual(self.proc._process_pulseintegrals()[0][1], -999.)
        self.assertAlmostEqual(self.proc._process_pulseintegrals()[0][3], 3.98951741969)
        self.proc.limit = None

    def create_tempfile_from_testdata(self):
        tmp_path = self.create_tempfile_path()
        data_path = self.get_testdata_path()
        shutil.copyfile(data_path, tmp_path)
        return tmp_path

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path

    def get_testdata_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_DATA_FILE)


class ProcessIndexedEventsTests(ProcessEventsTests):
    def setUp(self):
        warnings.filterwarnings('ignore')
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.open_file(self.data_path, 'a')
        self.proc = process_events.ProcessIndexedEvents(self.data, DATA_GROUP, [0, 10], progress=False)

    def test_process_traces(self):
        timings = self.proc.process_traces()
        self.assertEqual(timings[1][0], 162.5)
        self.assertEqual(timings[1][1], -999)

    def test_get_traces_for_indexed_event_index(self):
        self.assertEqual(self.proc.get_traces_for_indexed_event_index(0)[12][3], 1334)


class ProcessEventsWithLINTTests(ProcessEventsTests):
    def setUp(self):
        warnings.filterwarnings('ignore')
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.open_file(self.data_path, 'a')
        self.proc = process_events.ProcessEventsWithLINT(self.data, DATA_GROUP, progress=False)

    def test__reconstruct_time_from_traces(self):
        event = self.proc.source[10]
        times = self.proc._reconstruct_time_from_traces(event)
        self.assertAlmostEqual(times[0], 160.685483871)
        self.assertEqual(times[2], -999)

    def test__reconstruct_time_from_trace(self):
        trace = [200, 220]
        self.assertEqual(self.proc._reconstruct_time_from_trace(trace, 180), 0)
        self.assertEqual(self.proc._reconstruct_time_from_trace(trace, 190), 0.5)
        self.assertEqual(self.proc._reconstruct_time_from_trace(trace, 200), 1)
        self.assertEqual(self.proc._reconstruct_time_from_trace(trace, 210), -999)


class ProcessEventsWithTriggerOffsetTests(ProcessEventsTests):
    def setUp(self):
        warnings.filterwarnings('ignore')
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.open_file(self.data_path, 'a')
        self.proc = process_events.ProcessEventsWithTriggerOffset(self.data, DATA_GROUP, progress=False)

    def test__reconstruct_time_from_traces(self):
        event = self.proc.source[10]
        times = self.proc._reconstruct_time_from_traces(event)
        self.assertEqual(times[0], 162.5)
        self.assertEqual(times[2], -999)
        self.assertEqual(times[4], 165)

    def test__reconstruct_time_from_traces_with_external(self):
        self.proc.trigger = [0, 0, 0, 1]
        event = self.proc.source[10]
        times = self.proc._reconstruct_time_from_traces(event)
        self.assertEqual(times[0], 162.5)
        self.assertEqual(times[2], -999)
        self.assertEqual(times[4], -999)

    def test__first_above_thresholds(self):
        # 2 detectors
        self.assertEqual(self.proc._first_above_thresholds((x for x in [200, 200, 900]), [300, 400], 900), [2, 2, -999])
        self.assertEqual(self.proc._first_above_thresholds((x for x in [200, 200, 400]), [300, 400], 400), [2, 2, -999])
        self.assertEqual(self.proc._first_above_thresholds((x for x in [200, 350, 450, 550]), [300, 400], 550), [1, 2, -999])
        # 4 detectors
        self.assertEqual(self.proc._first_above_thresholds((x for x in [200, 200, 900]), [300, 400, 500], 900), [2, 2, 2])
        self.assertEqual(self.proc._first_above_thresholds((x for x in [200, 200, 400]), [300, 400, 500], 400), [2, 2, -999])
        self.assertEqual(self.proc._first_above_thresholds((x for x in [200, 350, 450, 550]), [300, 400, 500], 550), [1, 2, 3])
        # No signal
        self.assertEqual(self.proc._first_above_thresholds((x for x in [200, 250, 200, 2000]), [300, 400, 500], 250), [-999, -999, -999])

    def test__first_value_above_threshold(self):
        trace = [200, 200, 300, 200]
        self.assertEqual(self.proc._first_value_above_threshold(trace, 200), (0, 200))
        self.assertEqual(self.proc._first_value_above_threshold(trace, 250), (2, 300))
        self.assertEqual(self.proc._first_value_above_threshold(trace, 250, 4), (6, 300))
        self.assertEqual(self.proc._first_value_above_threshold(trace, 500), (-999, 0))

    def test__reconstruct_trigger(self):
        self.proc.trigger = (0, 0, False, 0)
        low_idx = [-999, -999, -999, -999]
        high_idx = [-999, -999, -999, -999]
        result = -999
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        self.proc.trigger = (0, 0, True, 0)
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)

        # Standard two detector trigger
        self.proc.trigger = (2, 0, False, 0)
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        high_idx = [-999, -999, 10, -999]
        low_idx = [-999, -999, 3, -999]
        result = -999
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        low_idx = [-999, 0, 3, 2]
        result = 2
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        low_idx = [0, 2, 4, -999]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        low_idx = [0, 2, 3, -999]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        low_idx = [0, 2, -999, -999]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        low_idx = [-999, -999, 3, 6]
        result = 6
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)

        # Standard four detector trigger
        self.proc.trigger = (3, 2, True, 0)
        low_idx = [-999, -999, -999, -999]
        high_idx = [-999, -999, -999, -999]
        result = -999
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        # Trigger on low
        low_idx = [7, 4, 1, -999]
        high_idx = [-999, -999, -999, -999]
        result = 7
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        high_idx = [8, 5, -999, -999]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        high_idx = [8, 9, 2, -999]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        # Trigger on high
        high_idx = [-999, 5, 2, -999]
        result = 5
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)

        # Other triggers
        self.proc.trigger = (1, 2, False, 0)
        low_idx = [1, 3, 5, 7]
        high_idx = [2, 4, -999, -999]
        result = 5
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        self.proc.trigger = (3, 0, False, 0)
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        self.proc.trigger = (0, 2, False, 0)
        result = 4
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        self.proc.trigger = (0, 4, False, 0)
        result = -999
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)
        self.proc.trigger = (1, 3, False, 0)
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx), result)


class ProcessEventsFromSourceTests(ProcessEventsTests):
    def setUp(self):
        warnings.filterwarnings('ignore')
        self.source_path = self.create_tempfile_from_testdata()
        self.source_data = tables.open_file(self.source_path, 'r')
        self.dest_path = self.create_tempfile_path()
        self.dest_data = tables.open_file(self.dest_path, 'a')
        self.proc = process_events.ProcessEventsFromSource(
            self.source_data, self.dest_data, DATA_GROUP, DATA_GROUP)

    def tearDown(self):
        warnings.resetwarnings()
        self.source_data.close()
        os.remove(self.source_path)
        self.dest_data.close()
        os.remove(self.dest_path)

    def test_process_and_store_results(self):
        self.proc.process_and_store_results()


class ProcessEventsFromSourceWithTriggerOffsetTests(ProcessEventsFromSourceTests,
                                                    ProcessEventsWithTriggerOffsetTests):
    def setUp(self):
        warnings.filterwarnings('ignore')
        self.source_path = self.create_tempfile_from_testdata()
        self.source_data = tables.open_file(self.source_path, 'r')
        self.dest_path = self.create_tempfile_path()
        self.dest_data = tables.open_file(self.dest_path, 'a')
        self.proc = process_events.ProcessEventsFromSourceWithTriggerOffset(
            self.source_data, self.dest_data, DATA_GROUP, DATA_GROUP)


class ProcessEventsFromSourceWithTriggerOffsetStationTests(ProcessEventsFromSourceTests,
                                                           ProcessEventsWithTriggerOffsetTests):
    def setUp(self):
        warnings.filterwarnings('ignore')
        self.source_path = self.create_tempfile_from_testdata()
        self.source_data = tables.open_file(self.source_path, 'r')
        self.dest_path = self.create_tempfile_path()
        self.dest_data = tables.open_file(self.dest_path, 'a')
        self.proc = process_events.ProcessEventsFromSourceWithTriggerOffset(
            self.source_data, self.dest_data, DATA_GROUP, DATA_GROUP,
            station=501)

    def test__reconstruct_time_from_traces_with_external(self):
        mock_trigger = Mock()
        mock_trigger.return_value = ([(process_events.ADC_LOW_THRESHOLD,
                                       process_events.ADC_HIGH_THRESHOLD)] * 4,
                                     [0, 0, 0, 1])
        self.proc.station.trigger = mock_trigger

        event = self.proc.source[10]
        times = self.proc._reconstruct_time_from_traces(event)
        self.assertEqual(times[0], 162.5)
        self.assertEqual(times[2], -999)
        self.assertEqual(times[4], -999)


if __name__ == '__main__':
    unittest.main()
