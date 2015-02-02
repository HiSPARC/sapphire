import unittest
import tempfile
import os
import shutil
import warnings

import tables

from sapphire.analysis import process_events


TEST_DATA_FILE = 'PE-testdata.h5'
DATA_GROUP = '/s501'


@unittest.skipIf(not os.path.exists(os.path.join(os.path.dirname(__file__), TEST_DATA_FILE)),
                 'Missing test datafile.')
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

    def test__reconstruct_time_from_traces(self):
        event = self.proc.source[10]
        times = self.proc._reconstruct_time_from_traces(event)
        self.assertEqual(times[0], 162.5)
        self.assertEqual(times[2], -999)

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
        # 2 detectors
        high_idx = []
        low_idx = [-999, 3]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx, n_detectors=2), -999)
        low_idx = [0, 3]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx, n_detectors=2), 3)

        # 4 detectors
        high_idx = [-999, -999]
        low_idx = [-999, 3]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx, n_detectors=4), -999)
        low_idx = [0, 3, 2]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx, n_detectors=4), 3)
        high_idx = [0, 3]
        low_idx = [0, 2, 4]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx, n_detectors=4), 3)
        high_idx = [0, 5]
        low_idx = [0, 2, 3]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx, n_detectors=4), 3)
        high_idx = [0, 5]
        low_idx = [0, 2]
        self.assertEqual(self.proc._reconstruct_trigger(low_idx, high_idx, n_detectors=4), 5)


if __name__ == '__main__':
    unittest.main()
