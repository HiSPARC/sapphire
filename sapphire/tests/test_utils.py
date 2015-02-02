import unittest
import types
from StringIO import StringIO

from numpy import pi, random
import progressbar

from sapphire import utils


class PbarTests(unittest.TestCase):

    def setUp(self):
        self.iterable = range(10)
        self.output = StringIO()

    def test_pbar_iterable(self):
        pb = utils.pbar(self.iterable, fd=self.output)
        self.assertIsInstance(pb, progressbar.ProgressBar)
        self.assertEqual(list(pb), self.iterable)

    def test_pbar_generator(self):
        """Return original generator, not a progressbar"""

        generator = (x for x in self.iterable)
        pb = utils.pbar(generator)
        self.assertIsInstance(pb, types.GeneratorType)
        self.assertEqual(list(pb), self.iterable)

    def test_pbar_generator_known_length(self):
        """Return progressbar for generator with known length"""

        generator = (y for y in self.iterable)
        pb = utils.pbar(generator, length=len(self.iterable), fd=self.output)
        self.assertIsInstance(pb, progressbar.ProgressBar)
        self.assertEqual(list(pb), self.iterable)

    def test_pbar_generator_wrong_length(self):
        """Raise exception for generator with wrong length"""

        generator = (y for y in self.iterable)
        pb = utils.pbar(generator, length=len(self.iterable) - 5, fd=self.output)
        self.assertRaises(ValueError, list, pb)

    def test_pbar_hide_output(self):
        """Empty output when not showing progressbar"""

        pb = utils.pbar(self.iterable, show=False, fd=self.output)
        self.assertEqual(list(pb), self.iterable)
        self.assertEqual(self.output.getvalue(), '')

        pb = utils.pbar(self.iterable, show=True, fd=self.output)
        self.assertEqual(list(pb), self.iterable)
        self.assertNotEqual(self.output.getvalue(), '')


class InBaseTests(unittest.TestCase):

    def test_ceil(self):
        self.assertEqual(utils.ceil_in_base(2.4, 2.5), 2.5)

    def test_floor(self):
        self.assertEqual(utils.floor_in_base(2.4, 2.5), 0)


class AngleBetweenTests(unittest.TestCase):

    def test_zeniths(self):
        for zenith in random.uniform(0, pi / 2, 10):
            self.assertAlmostEqual(utils.angle_between(zenith, 0, 0, 0), zenith)
            self.assertAlmostEqual(utils.angle_between(0, 0, zenith, 0), zenith)

    def test_azimuths(self):
        # Set both zeniths to pi/2 to give azimuth full effect.
        z = pi / 2
        for azimuth in random.uniform(-pi, pi, 20):
            self.assertAlmostEqual(utils.angle_between(z, azimuth, z, 0), abs(azimuth))
            self.assertAlmostEqual(utils.angle_between(z, 0, z, azimuth), abs(azimuth))

    def test_no_zenith(self):
        for azimuth in random.uniform(-pi, pi, 20):
            self.assertAlmostEqual(utils.angle_between(0, azimuth, 0, 0), 0)
            self.assertAlmostEqual(utils.angle_between(0, 0, 0, azimuth), 0)


if __name__ == '__main__':
    unittest.main()
