import unittest
import random
from math import pi

from sapphire.transformations import angles


class DegreeRadianHourTests(unittest.TestCase):

    def setUp(self):
        x = random.random()
        # (degrees, radians, hours)
        self.combinations = ((0., 0., 0.),
                             (15., pi / 12., 1.),
                             (90., pi / 2., 6.),
                             (180., pi, 12.),
                             (360. * x, 2 * pi * x, 24. * x))

    def test_hours_to_degrees(self):
        for degree, _, hour in self.combinations:
            self.assertAlmostEqual(angles.hours_to_degrees(hour), degree, 10)
            self.assertAlmostEqual(angles.degrees_to_hours(degree), hour, 10)
            self.assertAlmostEqual(angles.hours_to_degrees(angles.degrees_to_hours(degree)), degree, 10)
            self.assertAlmostEqual(angles.degrees_to_hours(angles.hours_to_degrees(hour)), hour, 10)

    def test_hours_to_radians(self):
        for _, radian, hour in self.combinations:
            self.assertAlmostEqual(angles.hours_to_radians(hour), radian)
            self.assertAlmostEqual(angles.radians_to_hours(radian), hour)
            self.assertAlmostEqual(angles.hours_to_radians(angles.radians_to_hours(radian)), radian)
            self.assertAlmostEqual(angles.radians_to_hours(angles.hours_to_radians(hour)), hour)


if __name__ == '__main__':
    unittest.main()
