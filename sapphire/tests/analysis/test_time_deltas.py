import unittest
import tempfile
import os
import shutil

import tables
from mock import patch, sentinel, Mock

from sapphire.analysis import time_deltas


TEST_DATA_FILE = 'test_data/esd_coincidences.h5'


class ProcessTimeDeltasTests(unittest.TestCase):
    def setUp(self):
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.open_file(self.data_path, 'a')
        self.td = time_deltas.ProcessTimeDeltas(self.data, progress=False)

    def tearDown(self):
        self.data.close()
        os.remove(self.data_path)

    def test_init(self):
        self.assertEqual(self.td.progress, False)
        self.assertEqual(self.td.data, self.data)

    def test_find_station_pairs(self):
        self.td.find_station_pairs()
        self.assertEqual(self.td.pairs, {(501, 502)})

    @patch.object(time_deltas, 'Station')
    def test_get_detector_offsets(self, mock_station):
        mock_offsets = Mock()
        mock_station.return_value = mock_offsets

        self.td.pairs = {(sentinel.station1, sentinel.station2),
                         (sentinel.station1, sentinel.station3)}
        self.td.get_detector_offsets()

        self.assertEqual(self.td.detector_timing_offsets,
                         {sentinel.station1: mock_offsets.detector_timing_offset,
                          sentinel.station2: mock_offsets.detector_timing_offset,
                          sentinel.station3: mock_offsets.detector_timing_offset})

    def test_store_time_deltas(self):
        pair = (501, 502)
        node_path = '/time_deltas/station_%d/station_%d' % pair
        self.assertRaises(Exception, self.data.get_node, node_path, 'time_deltas')
        self.td.store_time_deltas([12345678987654321], [2.5], pair)
        stored_data = self.data.get_node(node_path, 'time_deltas')
        self.assertEqual(list(stored_data[0]),
                         [12345678987654321, 12345678, 987654321, 2.5])

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


if __name__ == '__main__':
    unittest.main()
