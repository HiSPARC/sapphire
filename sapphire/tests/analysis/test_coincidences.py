import unittest
import tempfile
import os
import shutil

from mock import sentinel, patch, Mock
import tables
from numpy import uint64

from sapphire.analysis import coincidences
from sapphire.tests.validate_results import validate_results


TEST_DATA_FILE = 'test_data/esd_coincidences.h5'


class CoincidencesTests(unittest.TestCase):

    def setUp(self):
        self.c = coincidences.Coincidences(sentinel.data, None,
                                           sentinel.station_groups,
                                           progress=False)

    @patch.object(coincidences.Coincidences, 'search_coincidences')
    @patch.object(coincidences.Coincidences, 'process_events')
    @patch.object(coincidences.Coincidences, 'store_coincidences')
    def test_search_and_store_coincidences(self, mock_store, mock_process,
                                           mock_search):
        self.c.search_and_store_coincidences()
        mock_search.assert_called_with(window=2000)
        mock_process.assert_called_with()
        mock_store.assert_called_with(None)
        self.c.search_and_store_coincidences(sentinel.window, sentinel.cluster)
        mock_search.assert_called_with(window=sentinel.window)
        mock_process.assert_called_with()
        mock_store.assert_called_with(sentinel.cluster)

    def test__retrieve_timestamps(self):
        station1 = Mock()
        station2 = Mock()
        # Station 2 timestamps are not already correctly sorted.
        station1.col.return_value = [uint64(1400000002000000600), uint64(1400000008000000050)]
        station2.col.return_value = [uint64(1400000002000000700), uint64(1400000008000000000)][::-1]
        stations = [station1, station2]
        timestamps = self.c._retrieve_timestamps(stations)
        self.assertEqual(timestamps,
                         [(1400000002000000600, 0, 0), (1400000002000000700, 1, 1),
                          (1400000008000000000, 1, 0), (1400000008000000050, 0, 1)])
        # Shift both
        timestamps = self.c._retrieve_timestamps(stations, shifts=[-50, 10])
        self.assertEqual(timestamps,
                         [(1400000002000000550, 0, 0), (1400000002000000710, 1, 1),
                          (1400000008000000000, 0, 1), (1400000008000000010, 1, 0)])
        # Wrong value type shifts
        self.assertRaises(ValueError, self.c._retrieve_timestamps, stations, shifts=['', ''])
        self.assertRaises(ValueError, self.c._retrieve_timestamps, stations, shifts=['', 90])
        # Different length shifts
        timestamps = self.c._retrieve_timestamps(stations, shifts=[110])
        self.assertEqual(timestamps,
                         [(1400000002000000700, 1, 1), (1400000002000000710, 0, 0),
                          (1400000008000000000, 1, 0), (1400000008000000160, 0, 1)])
        timestamps = self.c._retrieve_timestamps(stations, shifts=[None, 60])
        self.assertEqual(timestamps,
                         [(1400000002000000600, 0, 0), (1400000002000000760, 1, 1),
                          (1400000008000000050, 0, 1), (1400000008000000060, 1, 0)])
        # Subnanosecond shifts
        timestamps = self.c._retrieve_timestamps(stations, shifts=[0.3, 5.9])
        self.assertEqual(timestamps,
                         [(1400000002000000600, 0, 0), (1400000002000000705, 1, 1),
                          (1400000008000000005, 1, 0), (1400000008000000050, 0, 1)])
        # Using limits
        timestamps = self.c._retrieve_timestamps(stations, limit=1)
        self.assertEqual(timestamps,
                         [(1400000002000000600, 0, 0), (1400000008000000000, 1, 0)])

    def test__do_search_coincidences(self):
        # [(timestamp, station_idx, event_idx), ..]
        timestamps = [(uint64(0), 0, 0), (uint64(0), 1, 0), (uint64(10), 1, 1),
                      (uint64(15), 2, 0), (uint64(100), 1, 2), (uint64(200), 2, 1),
                      (uint64(250), 0, 1), (uint64(251), 0, 2)]

        c = self.c._do_search_coincidences(timestamps, window=6)
        expected_coincidences = [[0, 1], [2, 3], [6, 7]]
        self.assertEqual(c, expected_coincidences)

        c = self.c._do_search_coincidences(timestamps, window=150)
        expected_coincidences = [[0, 1, 2, 3, 4], [4, 5], [5, 6, 7]]
        self.assertEqual(c, expected_coincidences)

        c = self.c._do_search_coincidences(timestamps, window=300)
        expected_coincidences = [[0, 1, 2, 3, 4, 5, 6, 7]]
        self.assertEqual(c, expected_coincidences)


class CoincidencesESDTests(CoincidencesTests):

    def setUp(self):
        self.c = coincidences.CoincidencesESD(sentinel.data, None,
                                              sentinel.station_groups,
                                              progress=False)

    @patch.object(coincidences.CoincidencesESD, 'search_coincidences')
    @patch.object(coincidences.CoincidencesESD, 'store_coincidences')
    def test_search_and_store_coincidences(self, mock_store, mock_search):
        self.c.search_and_store_coincidences()
        mock_search.assert_called_with(window=2000)
        mock_store.assert_called_with(cluster=None)
        self.c.search_and_store_coincidences(sentinel.window,
                                             sentinel.cluster)
        mock_search.assert_called_with(window=sentinel.window)
        mock_store.assert_called_with(cluster=sentinel.cluster)

    @patch.object(coincidences.CoincidencesESD, '_search_coincidences')
    def test_search_coincidences(self, mock__search):
        mock__search.return_value = (sentinel.c_index, sentinel.timestamps)
        self.c.search_coincidences()
        mock__search.assert_called_with(2000, None, None)
        self.assertEqual(self.c._src_timestamps, sentinel.timestamps)
        self.assertEqual(self.c._src_c_index, sentinel.c_index)

        self.c.search_coincidences(sentinel.window, sentinel.shifts,
                                   sentinel.limit)
        mock__search.assert_called_with(sentinel.window, sentinel.shifts,
                                        sentinel.limit)


class CoincidencesESDDataTests(unittest.TestCase):

    def setUp(self):
        self.data_path = self.create_tempfile_from_testdata()

    def tearDown(self):
        os.remove(self.data_path)

    def test_coincidencesesd_output(self):
        with tables.open_file(self.data_path, 'a') as data:
            c = coincidences.CoincidencesESD(data, '/coincidences',
                                             ['/station_501', '/station_502'],
                                             progress=False)
            c.search_and_store_coincidences()
        validate_results(self, self.get_testdata_path(), self.data_path)

    def create_tempfile_from_testdata(self):
        tmp_path = self.create_tempfile_path()
        data_path = self.get_testdata_path()
        shutil.copyfile(data_path, tmp_path)
        self.remove_existing_coincidences(tmp_path)
        return tmp_path

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path

    def get_testdata_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_DATA_FILE)

    def remove_existing_coincidences(self, path):
        with tables.open_file(path, 'a') as data:
            data.remove_node('/coincidences', recursive=True)


if __name__ == '__main__':
    unittest.main()
