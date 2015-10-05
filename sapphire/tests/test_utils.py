import unittest
import types
from StringIO import StringIO

from numpy import pi, random, exp, sqrt
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
        self.assertEqual(utils.ceil_in_base(0.1, 2.5), 2.5)

    def test_floor(self):
        self.assertEqual(utils.floor_in_base(2.4, 2.5), 0)
        self.assertEqual(utils.floor_in_base(0.1, 2.5), 0)

    def test_round(self):
        self.assertEqual(utils.round_in_base(2.4, 2.5), 2.5)
        self.assertEqual(utils.round_in_base(0.1, 2.5), 0)

    def test_zero_base(self):
        self.assertRaises(ZeroDivisionError, utils.ceil_in_base, 0.1, 0)
        self.assertRaises(ZeroDivisionError, utils.floor_in_base, 0.1, 0)
        self.assertRaises(ZeroDivisionError, utils.round_in_base, 0.1, 0)

    def test_integers(self):
        self.assertEqual(utils.ceil_in_base(3, 4), 4)
        self.assertEqual(utils.floor_in_base(3, 4), 0)
        self.assertEqual(utils.round_in_base(3, 4), 4)


class ActiveIndexTests(unittest.TestCase):

    def test_get_active_index(self):
        """Test if the bisection returns the correct index

        - If timestamp is before the first timestamp return index for
          first item
        - If timestamp is after last timestamp return index for last item
        - If timestamp is in the range return index of rightmost value
          equal or less than the timestamp

        """
        timestamps = [1., 2., 3., 4.]

        for idx, ts in [(0, 0.), (0, 1.), (0, 1.5), (1, 2.), (1, 2.1), (3, 4.),
                        (3, 5.)]:
            self.assertEqual(utils.get_active_index(timestamps, ts), idx)


class GaussTests(unittest.TestCase):

    """Test against explicit Gaussian"""

    def gaussian(self, x, N, mu, sigma):
        return N * exp(-(x - mu) ** 2. / (2. * sigma ** 2)) / (sigma * sqrt(2 * pi))

    def test_gauss(self):
        x, N, mu, sigma = (1., 1., 0., 1.)
        self.assertEqual(utils.gauss(x, N, mu, sigma), self.gaussian(x, N, mu, sigma))
        N = 2.
        self.assertEqual(utils.gauss(x, N, mu, sigma), self.gaussian(x, N, mu, sigma))
        sigma = 2.
        self.assertEqual(utils.gauss(x, N, mu, sigma), self.gaussian(x, N, mu, sigma))
        x = 1e5
        self.assertEqual(utils.gauss(x, N, mu, sigma), 0.)

    def test_gauss_array(self):
        """Test for arrays of random values"""

        n = 10000
        x, N, mu = random.uniform(-100, 100, size=(3, n))
        # sigma can not be 0
        sigma = random.uniform(1e-15, 100, size=n)
        value1 = utils.gauss(x, N, mu, sigma)
        value2 = self.gaussian(x, N, mu, sigma)
        self.assertTrue(all(abs(value1 - value2) < 1e-10))


class AngleBetweenTests(unittest.TestCase):

    """Check opening angle between two directions"""

    def test_zeniths(self):
        """One of the directions is the Zenith"""

        n = 10000
        zenith = random.uniform(0, pi / 2, n)
        azimuth1 = random.uniform(-pi, pi, n)
        azimuth2 = random.uniform(-pi, pi, n)
        angle = utils.angle_between(zenith, azimuth1, 0, azimuth2)
        self.assertTrue(all(abs(angle - zenith) < 1e-15))
        angle = utils.angle_between(0, azimuth1, zenith, azimuth2)
        self.assertTrue(all(abs(angle - zenith) < 1e-15))

    def test_azimuths(self):
        """Both directions at the horizon"""

        zenith = pi / 2
        azimuth = random.uniform(-pi, pi, 10000)
        angle = utils.angle_between(zenith, azimuth, zenith, 0)
        self.assertTrue(all(abs(angle - abs(azimuth)) < 1e-10))
        angle = utils.angle_between(zenith, 0, zenith, azimuth)
        self.assertTrue(all(abs(angle - abs(azimuth)) < 1e-10))

    def test_no_zenith(self):
        """Azimuths are irrelevant when from the Zenith"""

        azimuth1 = random.uniform(-pi, pi, 10000)
        azimuth2 = random.uniform(-pi, pi, 10000)
        angle = utils.angle_between(0, azimuth1, 0, azimuth2)
        self.assertTrue(all(angle == 0))

    def test_single_values(self):
        """Other tests use arrays, check if single values also work"""

        zenith = random.uniform(0, pi / 2)
        azimuth = random.uniform(-pi, pi)
        angle = utils.angle_between(zenith, azimuth, zenith, azimuth)
        self.assertTrue(angle == 0)


class DistanceBetweenTests(unittest.TestCase):

    """Check distance between two (x, y) cartesian coordinates"""

    def test_distances(self):
        """Check if distances are correctly calculated"""

        combinations = [((0, 0, 1.6, 0), 1.6),
                        ((-1, 0, 1, 0), 2),
                        ((-1, 0, -1, 0), 0),
                        ((random.uniform(1e-15, 100),) * 4, 0),
                        ((-10, -10, 5, 5), sqrt(450))]
        for coordinates, distance in combinations:
            self.assertEqual(utils.distance_between(*coordinates), distance)
            # same result if the coordinates and x, y are swapped
            self.assertEqual(utils.distance_between(*coordinates[::-1]), distance)


class WhichTests(unittest.TestCase):

    """Check if which works"""

    def test_which(self):
        """Check existence of common command"""

        utils.which('ls')

    def test_non_existent_program(self):
        """Check for error for non-existent program"""

        self.assertRaises(Exception, utils.which,
                          'a_very_unlikely_program_name_to_exist_cosmic_ray')


if __name__ == '__main__':
    unittest.main()
