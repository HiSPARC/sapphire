import unittest
import random
from math import pi

from sapphire.transformations import celestial


class ZenithAzimuthHorizontalTests(unittest.TestCase):

    def setUp(self):
        self.zenith = (0., pi / 4., pi / 2.)
        self.altitude = (pi / 2., pi / 4., 0.)

        self.azimuth = (-pi / 2., 0., pi / 2.)  # -pi
        self.Azimuth = (-pi, pi / 2., 0.)  # -pi / 2.

    def test_zenithazimuth_to_horizontal(self):
        for zenith, altitude in zip(self.zenith, self.altitude):
            self.assertEqual(celestial.zenithazimuth_to_horizontal(zenith, 0)[0], altitude)
            self.assertEqual(celestial.horizontal_to_zenithazimuth(altitude, 0)[0], zenith)

        for azimuth, Azimuth in zip(self.azimuth, self.Azimuth):
            self.assertEqual(celestial.zenithazimuth_to_horizontal(0, azimuth)[1], Azimuth)
            self.assertEqual(celestial.horizontal_to_zenithazimuth(0, Azimuth)[1], azimuth)


if __name__ == '__main__':
    unittest.main()
