import unittest

from sapphire.analysis import core_reconstruction


class BaseAlgorithm(object):

    """Use this class to check the different algorithms

    They should give similar results and errors in some cases.

    """

    def call_reconstruct(self, t, x, y, z):
        return self.algorithm.reconstruct_common(t, x, y, z)

    def test_stations_square(self):
        """Four detection points in a square shape."""

        # Same density
        p = (1., 1., 1., 1.)
        x = (0., 0., 10., 10.)
        y = (0., 10., 10., 0.)
        z = (0., 0., 0., 0.)
        result = self.call_reconstruct(p, x, y, z)
        self.assertAlmostEqual(result[0], 5.)
        self.assertAlmostEqual(result[1], 5.)


class CenterMassAlgorithmTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = core_reconstruction.CenterMassAlgorithm()


class AverageIntersectionAlgorithmTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = core_reconstruction.AverageIntersectionAlgorithm()


class EllipsLdfAlgorithmTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = core_reconstruction.EllipsLdfAlgorithm()


if __name__ == '__main__':
    unittest.main()
