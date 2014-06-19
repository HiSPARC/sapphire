from mock import sentinel, Mock, patch, call
import unittest

from numpy import nan, isnan, pi, degrees

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
        """Multiple detections at same point make reconstruction impossible."""

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

    def test_shower_at_angle(self):
        """Simple shower from a specific angle."""

        t = (0., 9., 9.)
        x = (0., 10., 0.)
        y = (0., 0., 10.)
        z = (0., 0., 0.)
        theta, phi = self.call_reconstruct(t, x, y, z)
        self.assertAlmostEqual(degrees(theta), 22.4476, 4)
        self.assertAlmostEqual(degrees(phi), -135.0000, 4)


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
