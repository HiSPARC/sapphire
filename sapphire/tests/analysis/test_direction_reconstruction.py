from mock import sentinel, Mock, patch, call
import unittest

from numpy import nan, isnan, pi, degrees

from sapphire.analysis import direction_reconstruction


class BaseAlgorithm(object):

    """Use this class to check the different algorithms

    They should give similar results and errors in some cases.

    """

    def call_reconstruct(self, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0, z1, z2):
        return self.algorithm.reconstruct_common(t0, t1, t2, x0, x1, x2,
                                                 y0, y1, y2, z0, z1, z2)

    def test_stations_in_line(self):
        """Three detection points on a line does not provide a solution."""

        # On a line in x
        t0, t1, t2 = (0., 2., 3.)
        x0, x1, x2 = (0., 0., 0.)  # same x
        y0, y1, y2 = (0., 5., 10.)
        z0, z1, z2 = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t0, t1, t2, x0, x1, x2, y0, y1, y2,
                                       z0, z1, z2)
        self.assertTrue(isnan(result).all())

        # Diagonal line
        t0, t1, t2 = (0., 2., 3.)
        x0, x1, x2 = (0., 5., 10.)
        y0, y1, y2 = (0., 5., 10.)
        z0, z1, z2 = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t0, t1, t2, x0, x1, x2, y0, y1, y2,
                                       z0, z1, z2)
        self.assertTrue(isnan(result).all())

    def test_same_stations(self):
        """Multiple detections at same point make reconstruction impossible."""

        # Two at same location
        t0, t1, t2 = (0., 2., 3.)
        x0, x1, x2 = (0., 0., 1.)
        y0, y1, y2 = (5., 5., 6.)
        z0, z1, z2 = (0., 0., 1.)
        result = self.call_reconstruct(t0, t1, t2, x0, x1, x2, y0, y1, y2,
                                       z0, z1, z2)
        self.assertTrue(isnan(result).all())

        # Three at same location
        t0, t1, t2 = (0., 2., 3.)
        x0, x1, x2 = (0., 0., 0.)  # same x
        y0, y1, y2 = (5., 5., 5.)  # same y
        z0, z1, z2 = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t0, t1, t2, x0, x1, x2, y0, y1, y2,
                                       z0, z1, z2)
        self.assertTrue(isnan(result).all())

    def test_shower_from_above(self):
        """Simple shower from zenith, azimuth can be any allowed value."""

        t0, t1, t2 = (0., 0., 0.)  # same t
        x0, x1, x2 = (0., 10., 0.)
        y0, y1, y2 = (0., 0., 10.)
        z0, z1, z2 = (0., 0., 0.)  # same z
        theta, phi = self.call_reconstruct(t0, t1, t2, x0, x1, x2, y0, y1, y2,
                                           z0, z1, z2)
        self.assertEqual(theta, 0)
        # azimuth can be any value between -pi and pi
        self.assertTrue(-pi <= phi <= pi)

    def test_shower_at_angle(self):
        """Simple shower from zenith, azimuth can be any allowed value."""

        t0, t1, t2 = (0., 9., 9.)  # same t
        x0, x1, x2 = (0., 10., 0.)
        y0, y1, y2 = (0., 0., 10.)
        z0, z1, z2 = (0., 0., 0.)  # same z
        theta, phi = self.call_reconstruct(t0, t1, t2, x0, x1, x2, y0, y1, y2,
                                           z0, z1, z2)
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
