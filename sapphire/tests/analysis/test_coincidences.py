from mock import sentinel, Mock, patch, call
import unittest

from sapphire.analysis import coincidences


class CoincidencesTests(unittest.TestCase):

    def setUp(self):
        self.coincidences = coincidences.Coincidences(sentinel.data,
                                                      None,
                                                      sentinel.station_groups,
                                                      progress=False)

    def test__do_search_coincidences(self):
        timestamps = [(0, 0, 0), (0, 1, 0), (10, 1, 1), (15, 2, 0),
                      (100, 1, 2), (200, 2, 1), (250, 0, 1), (251, 0, 2)]

        c = self.coincidences._do_search_coincidences(timestamps, window=6)
        expected_coincidences = [[0, 1], [2, 3], [6, 7]]
        self.assertEqual(c, expected_coincidences)

        c = self.coincidences._do_search_coincidences(timestamps, window=150)
        expected_coincidences = [[0, 1, 2, 3, 4], [4, 5], [5, 6, 7]]
        self.assertEqual(c, expected_coincidences)

        c = self.coincidences._do_search_coincidences(timestamps, window=300)
        expected_coincidences = [[0, 1, 2, 3, 4, 5, 6, 7]]
        self.assertEqual(c, expected_coincidences)


if __name__ == '__main__':
    unittest.main()
