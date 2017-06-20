import calendar
import datetime
import unittest

from math import pi


import numpy as np

from sapphire.transformations import base, celestial, clock


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


class EquatorialTests(unittest.TestCase):

    """Accuracy tests for Celestial coordinate transformations.

    Use references as tests also used by astropy.

    Source:
    https://github.com/astropy/astropy/blob/master/astropy/coordinates/tests/accuracy/test_altaz_icrs.py

    Licensed under a 3-clause BSD style license - see:
    https://github.com/astropy/astropy/blob/master/licenses/LICENSE.rst

    """

    def test_against_hor2eq(self):
        """Check for consistent results with an IDL hor2eq example.

        See EXAMPLE input and output here:
        http://idlastro.gsfc.nasa.gov/ftp/pro/astro/hor2eq.pro

        Observatory position for ``kpno`` from here:
        http://idlastro.gsfc.nasa.gov/ftp/pro/astro/observatory.pro

        """
        # IDL hor2eq
        ra_expected = np.radians(3.30875)
        dec_expected = np.radians(15.183416666666666)

        # Astropy 1.0rc1
        ra_astropy = np.radians(3.3094193224314625)
        dec_astropy = np.radians(15.183757021354532)

        # KPNO observatory
        longitude = -111.6
        latitude = 31.9633

        # Observation time
        jd = 2466879.7083333

        # Altitude Azimuth
        altitude = (37, 54, 41)
        azi = (264, 55, 6)

        # lst = clock.gmst_to_lst(clock.juliandate_to_gmst(jd), longitude)
        # Matches  LAST = +03 53 53.6  in the hor2eq.pro

        # SAPPHiRE
        utc = calendar.timegm(clock.juliandate_to_utc(jd).utctimetuple())
        gps = clock.utc_to_gps(utc)
        zenith, azimuth = celestial.horizontal_to_zenithazimuth(
            np.radians(base.sexagesimal_to_decimal(*altitude)),
            np.radians(base.sexagesimal_to_decimal(*azi)))
        ra, dec = celestial.zenithazimuth_to_equatorial(latitude, longitude,
                                                        gps, zenith, azimuth)

        # Test eq_to_zenaz merely against IDL
        zencalc, azcalc = celestial.equatorial_to_zenithazimuth(latitude,
        longitude, gps, ra_expected, dec_expected)

        self.assertAlmostEqual(ra, ra_expected, 1)
        self.assertAlmostEqual(ra, ra_astropy, 1)
        self.assertAlmostEqual(dec, dec_expected, 2)
        self.assertAlmostEqual(dec, dec_astropy, 2)

        self.assertAlmostEqual(zencalc, zenith, 1)
        self.assertAlmostEqual(azcalc, azimuth, 2)

    def test_against_pyephem(self):
        """Check for consistent results with one PyEphem example.

        PyEphem: http://rhodesmill.org/pyephem/

        See example input and output here:
        https://gist.github.com/zonca/1672906
        https://github.com/phn/pytpm/issues/2#issuecomment-3698679

        """
        # PyEphem
        ra_expected = np.radians(196.497518)
        dec_expected = np.radians(-4.569323)

        # Astropy 1.0rc1
        ra_astropy = np.radians(196.49537283)
        dec_astropy = np.radians(-4.5606942763)

        # Data
        longitude = base.sexagesimal_to_decimal(-109, -24, -53.1)
        latitude = base.sexagesimal_to_decimal(33, 41, 46.0)
        utc = datetime.datetime(2011, 9, 18, 8, 50)
        altitude = np.radians(-60.7665)
        azi = np.radians(6.8927)

        # SAPPHiRE
        gps = clock.utc_to_gps(calendar.timegm(utc.utctimetuple()))
        zenith, azimuth = celestial.horizontal_to_zenithazimuth(altitude, azi)
        ra, dec = celestial.zenithazimuth_to_equatorial(latitude, longitude,
                                                        gps, zenith, azimuth)

        zencalc, azcalc = celestial.equatorial_to_zenithazimuth(latitude,
                                                                longitude, gps,
                                                                ra_expected,
                                                                dec_expected)


        self.assertAlmostEqual(ra, ra_expected, 2)
        self.assertAlmostEqual(ra, ra_astropy, 2)
        self.assertAlmostEqual(dec, dec_expected, 2)
        self.assertAlmostEqual(dec, dec_astropy, 2)

        self.assertAlmostEqual(zencalc, zenith, 2)
        self.assertAlmostEqual(azcalc, azimuth, 2)

    def test_against_jpl_horizons(self):
        """Check for consistent results with the JPL Horizons example.

        The input parameters and reference results are taken from this page:
        (from the first row of the Results table at the bottom of that page)
        http://ssd.jpl.nasa.gov/?horizons_tutorial

        """
        # NASA JPL
        ra_expected = np.radians(291.229208333)
        dec_expected = np.radians(-40.9413611111)

        # Astropy 1.0rc1
        ra_astropy = np.radians(291.229161499)
        dec_astropy = np.radians(-40.9413052259)

        # Data
        # Kitt Peak
        longitude = 248.405300
        latitude = 31.9585
        utc = datetime.datetime(1998, 7, 28, 3, 0)
        altitude = np.radians(2.6223)
        azi = np.radians(143.2970)

        # SAPPHiRE
        gps = clock.utc_to_gps(calendar.timegm(utc.utctimetuple()))
        zenith, azimuth = celestial.horizontal_to_zenithazimuth(altitude, azi)
        ra, dec = celestial.zenithazimuth_to_equatorial(latitude, longitude,
                                                        gps, zenith, azimuth)

        zencalc, azcalc = celestial.equatorial_to_zenithazimuth(latitude,
                                                                longitude, gps,
                                                                ra_expected,
                                                                dec_expected)

        self.assertAlmostEqual(ra, ra_expected, 3)
        self.assertAlmostEqual(ra, ra_astropy, 3)
        self.assertAlmostEqual(dec, dec_expected, 2)
        self.assertAlmostEqual(dec, dec_astropy, 2)

        self.assertAlmostEqual(zencalc, zenith, 2)
        self.assertAlmostEqual(azcalc, azimuth, 2)

class AstropyEquatorialTests(unittest.TestCase):
    """
    This tests the 4 new astropy functions. They should be very close to
    Pyephem results and in this test they are compared to 10 different
    coordinates from astropy.
    """
    def test_pyephem(self):
        from sapphire.transformations import celestial
        import numpy as np
        # This is the transform inputs
        eq = [(-39.34633914878846, -112.2277168069694, 1295503840,
               3.8662384455822716, -0.31222454326513827),
              (53.13143508448587,-49.24074935964933, 985619982,
               3.901575896592809, -0.3926720112815971),
              (48.02031016860923, -157.4023812557098, 1126251396,
               3.366278312183976, -1.3610394240813288)]
        # transform inputs converted to of zenaz inputs
        zeneq = [(-39.34633914878846, -112.2277168069694,
                  1295503840, 1.8830208700600348, -2.295442118787375),
                 (53.13143508448587, -49.24074935964933, 985619982,
                  1.9634683380764937, -2.3307795697979126),
                 (48.02031016860923, -157.4023812557098, 1126251396,
                  2.9318357508762256, -1.7954819853890793)]

        # result of pyephem hor->eq/zenaz-> eq
        efemeq = [(5.620508199785029, -0.3651173667585858),
                 (5.244630787139936, -0.7866376569183651),
                 (2.276751381056623, -1.0406498066785745)]
        # result of pyephem eq->hor
        altaz = [(2.175107479095459, -0.19537943601608276),
                 (5.25273323059082, -0.8308737874031067),
                 (3.4536221027374268, -0.894329845905304)]
        # result converted to eq->zenaz
        zenaz = [(1.7661757628109793, -0.6043111523005624),
                 (2.4016701141980032, 2.6012484033836625),
                 (2.4651261727002005, -1.8828257759425302)]

        htoea_test = []
        etoha_test = []
        zenaztoeq_test = []
        eqtozenaz_test = []

        # Produce horizontal_to_equatorial_astropy results
        for i in eq:
            result = celestial.horizontal_to_equatorial_astropy(
                i[0], i[1], i[2], [(i[3], i[4])])
            htoea_test.extend(result)

        # Produce equatorial_to_horizontal_astropy results
        for i in eq:
            result = celestial.equatorial_to_horizontal_astropy(i[0], i[1],
                                                                i[2],
                                                                [(i[3], i[4])])
            etoha_test.extend(result)

        # Produce zenithazimuth_to_equatorial_astropy results
        for i in zeneq:
            result = celestial.zenithazimuth_to_equatorial_astropy(
                i[0], i[1], i[2], [(i[3], i[4])])
            zenaztoeq_test.extend(result)

        # Produce equatorial_to_zenithazimuth_astropy results
        for i in eq:
            result = celestial.equatorial_to_zenithazimuth_astropy(i[0], i[1],
                                                                i[2],
                                                                [(i[3],
                                                                  i[4])])
            eqtozenaz_test.extend(result)

        # Check if all inputs are correct, cast to numpy array for certainty
        # Test horizontal_to_equatorial_astropy
        np.testing.assert_almost_equal(efemeq, htoea_test, 4)
        # Test zenithazimuth_to_equatorial_astropy
        np.testing.assert_almost_equal(efemeq, zenaztoeq_test, 4)
        # Test equatorial_to_horizontal_astropy
        np.testing.assert_almost_equal(altaz, etoha_test, 4)
        # Test equatorial_to_zenithazimuth_astropy
        np.testing.assert_almost_equal(zenaz, eqtozenaz_test, 4)





if __name__ == '__main__':
    unittest.main()
