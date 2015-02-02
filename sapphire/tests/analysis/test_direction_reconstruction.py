import unittest

from numpy import isnan, pi, sqrt, arcsin, arctan

from sapphire.analysis import direction_reconstruction


class BaseAlgorithm(object):

    """Use this class to check the different algorithms

    They should give similar results and errors in some cases.

    """

    def call_reconstruct(self, t, x, y, z):
        return self.algorithm.reconstruct_common(t, x, y, z)

    def test_stations_in_line(self):
        """Three detection points on a line do not provide a solution."""

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
        self.assertAlmostEqual(theta, 0., 4)
        # azimuth can be any value between -pi and pi
        self.assertTrue(-pi <= phi <= pi)

    def test_show_to_large_dt(self):
        """Time difference larger than expected by speed of light."""

        # TODO: Add better test with smaller tolerance

        x = (0., -5., 5.)
        y = (sqrt(100 - 25), 0., 0.)
        z = (0., 0., 0.)

        t = (35., 0., 0.)
        theta, phi = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(theta))

    def test_showers_at_various_angles(self):
        """Simple shower from specific zenith angles."""

        c = .3

        x = (0., -5., 5.)
        y = (sqrt(100 - 25), 0., 0.)
        z = (0., 0., 0.)

        # triangle height
        h = sqrt(100 - 25)

        times = (2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5)

        for time in times:
            for i in range(3):
                zenith = arcsin((time * c) / h)

                t = [0., 0., 0.]
                t[i] = time
                azimuths = [-pi / 2, pi / 6, pi * 5 / 6]
                theta, phi = self.call_reconstruct(t, x, y, z)
                self.assertAlmostEqual(phi, azimuths[i], 4)
                self.assertAlmostEqual(theta, zenith, 4)

                t = [time] * 3
                t[i] = 0.
                azimuths = [pi / 2, -pi * 5 / 6, -pi / 6]
                theta, phi = self.call_reconstruct(t, x, y, z)
                self.assertAlmostEqual(phi, azimuths[i], 4)
                self.assertAlmostEqual(theta, zenith, 4)


class AltitudeAlgorithm(object):

    """Use this class to check the altitude support

    They should give similar results and errors in some cases.

    """

    def call_reconstruct(self, t, x, y, z):
        return self.algorithm.reconstruct_common(t, x, y, z)

    def test_stations_altitude(self):
        """Simple shower on a non horizontal square."""

        x = (0., 10., 10.)
        y = (0, 0., 10.)
        z = (2., 0., -2.)

        zenith = arctan(4. / 10. / sqrt(2))

        t = [0., 0., 0.]
        azimuth = pi / 4.
        theta, phi = self.call_reconstruct(t, x, y, z)

        self.assertAlmostEqual(phi, azimuth, 5)
        self.assertAlmostEqual(theta, zenith, 5)


class MultiAlgorithm(object):

    """Use this class to check the different algorithms for more stations

    They should give similar results and errors in some cases.

    """

    def call_reconstruct(self, t, x, y, z):
        return self.algorithm.reconstruct_common(t, x, y, z)

    def test_diamond_stations(self):
        """Simple shower from specific zenith angles."""

        c = .3

        x = (0., -5., 5., 10.)
        y = (sqrt(100 - 25), 0., 0., sqrt(100 - 25))
        z = (0., 0., 0., 0.)

        # triangle height
        h = sqrt(100 - 25)

        times = (2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5)

        for time in times:
            zenith = arcsin((time * c) / h)

            t = [0., 0., 0., 0.]
            t[1] = time
            t[3] = -time
            azimuth = pi / 6
            theta, phi = self.call_reconstruct(t, x, y, z)
            self.assertAlmostEqual(phi, azimuth, 5)
            self.assertAlmostEqual(theta, zenith, 5)

    def test_square_stations(self):
        """Simple shower from specific zenith angles."""

        c = .3

        x = (0., 5., 5., 0.)
        y = (0, 0., 5., 5.)
        z = (0., 0., 0., 0.)

        # triangle height
        h = sqrt(50. / 4.)

        times = (2.5, 5., 7.5, 10.)

        for time in times:
            zenith = arcsin((time * c) / h)

            t = [0., 0., 0., 0.]
            t[0] = -time
            t[2] = time
            azimuth = - 3 * pi / 4
            theta, phi = self.call_reconstruct(t, x, y, z)
            self.assertAlmostEqual(phi, azimuth, 5)
            self.assertAlmostEqual(theta, zenith, 5)


class MultiAltitudeAlgorithm(MultiAlgorithm, AltitudeAlgorithm):

    """Check some algorithms for multiple stations at different altitudes.

    They should give similar results and errors in some cases.

    """

    def test_hexagon_altitude(self):
        """Simple shower on a non horizontal square."""

        x = (-5., 5., 10., 5., -5., -10.)
        y = (-5. * sqrt(3), -5. * sqrt(3), 0., 5. * sqrt(3), 5. * sqrt(3), 0.)
        z = (0., -3., -5., -3., 0., 4,)

        zenith = 0.38333
        azimuth = 0.00000

        t = [0., 0., 0., 0., 0., 0.]
        theta, phi = self.call_reconstruct(t, x, y, z)

        self.assertAlmostEqual(phi, azimuth, 4)
        self.assertAlmostEqual(theta, zenith, 4)


class DirectAlgorithmTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithm()


class DirectAlgorithmCartesian2DTest(unittest.TestCase, BaseAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian2D()


class DirectAlgorithmCartesian3DTest(unittest.TestCase, BaseAlgorithm,
                                     AltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian3D()


class FitAlgorithmTest(unittest.TestCase, BaseAlgorithm,
                       MultiAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.FitAlgorithm()


class RegressionAlgorithmTest(unittest.TestCase, BaseAlgorithm, MultiAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.RegressionAlgorithm()


class RegressionAlgorithm3DTest(unittest.TestCase, BaseAlgorithm,
                                MultiAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.RegressionAlgorithm3D()


if __name__ == '__main__':
    unittest.main()
