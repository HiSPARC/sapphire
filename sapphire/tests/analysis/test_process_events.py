import unittest
import tempfile
import os
import shutil
import tables
import sys

from mock import patch
from sapphire.analysis import process_events


TEST_DATA_FILE = 'PE-testdata.h5'
DATA_GROUP = '/s501'


@unittest.skipIf(not os.path.exists(TEST_DATA_FILE), 'Missing test datafile.')
class ProcessEventsTests(unittest.TestCase):
    def setUp(self):
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.openFile(self.data_path, 'a')
        self.proc = process_events.ProcessEvents(self.data, DATA_GROUP)

        # make progressbar(list) do nothing (i.e., return list)
        self.progressbar_patcher = patch('progressbar.ProgressBar')
        self.progressbar_mock = self.progressbar_patcher.start()
        self.progressbar_mock.return_value.side_effect = lambda x: x

    def tearDown(self):
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

    def test__first_above_threshold(self):
        trace = [0, 2, 4, 2, 0]
        self.assertEqual(self.proc._first_above_threshold(trace, 1), 1)
        self.assertEqual(self.proc._first_above_threshold(trace, 3), 2)
        self.assertEqual(self.proc._first_above_threshold(trace, 4), 2)
        self.assertEqual(self.proc._first_above_threshold(trace, 5), -999)

    def test__process_pulseintegrals(self):
        self.proc.limit = 1
        self.assertAlmostEqual(self.proc._process_pulseintegrals()[0][3], 4.37504635)
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
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.openFile(self.data_path, 'a')
        self.proc = process_events.ProcessIndexedEvents(self.data, DATA_GROUP, [0, 10])

    def test_process_traces(self):
        timings = self.proc.process_traces()
        self.assertEqual(timings[1][0], 162.5)
        self.assertEqual(timings[1][1], -999)

    def test_get_traces_for_indexed_event_index(self):
        self.assertEqual(self.proc.get_traces_for_indexed_event_index(0)[12][3], 1334)


class ProcessEventsWithLINTTests(ProcessEventsTests):
    def setUp(self):
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.openFile(self.data_path, 'a')
        self.proc = process_events.ProcessEventsWithLINT(self.data, DATA_GROUP)

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
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.openFile(self.data_path, 'a')
        self.proc = process_events.ProcessEventsWithTriggerOffset(self.data, DATA_GROUP)

    def test__reconstruct_trigger_time_from_traces(self):
        # 2 detectors
        traces = [[200, 300], [200, 300]]  # 2 low
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces), 1)
        traces = [[200, 300], [200, 200]]  # 1 low, no trigger
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces), -999)

        # 4 detectors
        traces = [[300, 200], [300, 200], [200, 300], [200, 200]]  # 3 low
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces), 1)
        traces = [[400, 200], [300, 400], [200, 200], [200, 200]]  # 2 high (first 2 low)
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces), 1)
        traces = [[400, 200], [400, 200], [200, 300], [200, 200]]  # first 2 high, then 3 low
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces), 0)
        traces = [[300, 200], [300, 400], [300, 400], [200, 200]]  # first 3 low, then 2 high
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces), 0)
        traces = [[200, 200], [200, 200], [200, 200], [200, 200]]  # no trigger
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces), -999)

        # 3 detectors?!
        traces = [[200, 200], [200, 200], [200, 200]]
        self.assertRaises(LookupError, self.proc._reconstruct_trigger_time_from_traces, traces)
        traces = [[400, 200], [400, 200], [200, 200]]  # process as if 4 detectors, 2 high
        self.assertEqual(self.proc._reconstruct_trigger_time_from_traces(traces, 4), 0)


if __name__ == '__main__':
    unittest.main()
