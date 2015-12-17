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

        self.assertAlmostEqual(ra, ra_expected, 1)
        self.assertAlmostEqual(ra, ra_astropy, 1)
        self.assertAlmostEqual(dec, dec_expected, 2)
        self.assertAlmostEqual(dec, dec_astropy, 2)

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

        self.assertAlmostEqual(ra, ra_expected, 2)
        self.assertAlmostEqual(ra, ra_astropy, 2)
        self.assertAlmostEqual(dec, dec_expected, 2)
        self.assertAlmostEqual(dec, dec_astropy, 2)

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

        self.assertAlmostEqual(ra, ra_expected, 3)
        self.assertAlmostEqual(ra, ra_astropy, 3)
        self.assertAlmostEqual(dec, dec_expected, 2)
        self.assertAlmostEqual(dec, dec_astropy, 2)


if __name__ == '__main__':
    unittest.main()
