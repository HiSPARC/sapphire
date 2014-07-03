from mock import sentinel, Mock, patch, call
import unittest

from numpy import nan, isnan, pi, degrees, sqrt

from sapphire.analysis import direction_reconstruction


class BaseAlgorithm(object):

    """Use this class to check the different algorithms

    They should give similar results and errors in some cases.

    """

    def call_reconstruct(self, t, x, y, z):
        return self.algorithm.reconstruct_common(t, x, y, z)

    def test_stations_in_line(self):
        """Three detection points on a line does not provide a solution."""

        # On a line in x
        t = (0., 2., 3.)
        x = (0., 0., 0.)  # same x
        y = (0., 5., 10.)
        z = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

        # Diagonal line
        t = (0., 2., 3.)
        x = (0., 5., 10.)
        y = (0., 5., 10.)
        z = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

    def test_same_stations(self):
        """Multiple detections at same point make reconstruction impossible.

        With different arrival time.

        """
        # Two at same location
        t = (0., 2., 3.)
        x = (0., 0., 1.)
        y = (5., 5., 6.)
        z = (0., 0., 1.)
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

        # Three at same location
        t = (0., 2., 3.)
        x = (0., 0., 0.)  # same x
        y = (5., 5., 5.)  # same y
        z = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

    def test_shower_from_above(self):
        """Simple shower from zenith, azimuth can be any allowed value."""

        t = (0., 0., 0.)  # same t
        x = (0., 10., 0.)
        y = (0., 0., 10.)
        z = (0., 0., 0.)  # same z
        theta, phi = self.call_reconstruct(t, x, y, z)
        self.assertEqual(theta, 0)
        # azimuth can be any value between -pi and pi
        self.assertTrue(-pi <= phi <= pi)

    def test_shower_at_azimuths(self):
        """Simple shower from specific azimuth angles."""

        x = (0., -5., 5.)
        y = (sqrt(100 - 25), 0., 0.)
        z = (0., 0., 0.)

        # combinations of times and azimuths
        combinations = (((10., 0., 0.), -90.000),
                        ((0., 10., 0.), 30.000),
                        ((0., 0., 10.), 150.000),
                        ((10., 10., 0.), -30.000),
                        ((10., 0., 10.), -150.000),
                        ((0., 10., 10.), 90.000))

        for t, azimuth in combinations:
            theta, phi = self.call_reconstruct(t, x, y, z)
            self.assertAlmostEqual(degrees(phi), azimuth, 3)
            self.assertAlmostEqual(degrees(theta), 20.268, 3)

    def test_shower_at_zeniths(self):
        """Simple shower from specific zenith angles."""

        x = (0., -5., 5.)
        y = (sqrt(100 - 25), 0., 0.)
        z = (0., 0., 0.)

        # combinations of times and zeniths
        combinations = (((2.5, 0., 0.), 4.968),
                        ((5., 0., 0.), 9.974),
                        ((7.5, 0., 0.), 15.059),
                        ((10., 0., 0.), 20.268),
                        ((12.5, 0., 0.), 25.659),
                        ((15., 0., 0.), 31.306),
                        ((17.5, 0., 0.), 37.317),
                        ((20., 0., 0.), 43.854),
                        ((22.5, 0., 0.), 51.208),
                        ((25., 0., 0.), 60.000),
                        ((27.5, 0., 0.), 72.294))

        for t, zenith in combinations:
            theta, phi = self.call_reconstruct(t, x, y, z)
            self.assertAlmostEqual(degrees(phi), -90.00, 3)
            self.assertAlmostEqual(degrees(theta), zenith, 3)

        t = (0., 0., 0.)
        theta, phi = self.call_reconstruct(t, x, y, z)
        self.assertTrue(-pi <= phi <= pi)
        self.assertAlmostEqual(degrees(theta), 0.000, 3)

        t = (35., 0., 0.)
        theta, phi = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(theta))


class DirectAlgorithmTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithm()


class DirectAlgorithmCartesian2DTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian2D()


class DirectAlgorithmCartesian3DTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian3D()


class FitAlgorithmTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.FitAlgorithm()


if __name__ == '__main__':
    unittest.main()
