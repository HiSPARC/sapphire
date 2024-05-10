import unittest

from sapphire.analysis import core_reconstruction


class BaseAlgorithm:
    """Use this class to check the different algorithms

    They should give similar results and errors in some cases.

    """

    def call_reconstruct(self, t, x, y, z):
        return self.algorithm.reconstruct_common(t, x, y, z)

    def test_stations_square(self):
        """Four detection points in a square shape."""

        # Same density
        p = (1.0, 1.0, 1.0, 1.0)
        x = (0.0, 0.0, 10.0, 10.0)
        y = (0.0, 10.0, 10.0, 0.0)
        z = (0.0, 0.0, 0.0, 0.0)
        result = self.call_reconstruct(p, x, y, z)
        self.assertAlmostEqual(result[0], 5.0)
        self.assertAlmostEqual(result[1], 5.0)


class CenterMassAlgorithmTest(unittest.TestCase, BaseAlgorithm):
    def setUp(self):
        self.algorithm = core_reconstruction.CenterMassAlgorithm()


class AverageIntersectionAlgorithmTest(unittest.TestCase, BaseAlgorithm):
    def setUp(self):
        self.algorithm = core_reconstruction.AverageIntersectionAlgorithm()


class EllipsLdfAlgorithmTest(unittest.TestCase, BaseAlgorithm):
    def setUp(self):
        self.algorithm = core_reconstruction.EllipsLdfAlgorithm()
