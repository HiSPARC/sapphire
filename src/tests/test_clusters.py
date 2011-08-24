import unittest
from math import pi, atan

from mock import Mock

import clusters

class DetectorTests(unittest.TestCase):
    def setUp(self):
        self.mock_station = Mock()
        self.detector_1 = clusters.Detector(self.mock_station, 1, 0, 'LR')
        self.detector_2 = clusters.Detector(self.mock_station, -1, 2, 'UD')

    def test_detector_size(self):
        self.assertEqual(self.detector_1.detector_size, (.5, 1.))
        self.assertEqual(self.detector_2.detector_size, (.5, 1.))

    def test_attributes(self):
        self.assertIs(self.detector_1.station, self.mock_station)
        self.assertEqual(self.detector_1.x, 1)
        self.assertEqual(self.detector_1.y, 0)
        self.assertEqual(self.detector_1.orientation, 'LR')

    def test_LR_get_corners(self):
        corners = self.detector_1.get_corners(0.25, 3, 0)
        self.assertEqual(corners, [(.75, 2.75), (1.75, 2.75), (1.75, 3.25), (.75, 3.25)])

    def test_LR_get_corners_rotated(self):
        corners = self.detector_1.get_corners(0, 0, pi / 2)
        expected_corners = [(.25, .5), (.25, 1.5), (-.25, 1.5), (-.25, .5)]
        for (x, y), (expected_x, expected_y) in zip(corners, expected_corners):
            self.assertAlmostEqual(x, expected_x)
            self.assertAlmostEqual(y, expected_y)

    def test_UD_get_corners(self):
        corners = self.detector_2.get_corners(0.25, 3, 0)
        self.assertEqual(corners, [(-1, 4.5), (-.5, 4.5), (-.5, 5.5), (-1, 5.5)])


if __name__ == '__main__':
    unittest.main()
