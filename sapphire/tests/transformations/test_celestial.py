import calendar
import datetime
import unittest
import random as r
from math import pi
import matplotlib.pyplot as plt

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


def test_transformations():
    # This tuple list contains all the used testcases for the test.
    # The meanings of the tuple entries are:
    # 0: latitude (deg) 1: longitude (deg) 2: utc_timestamp 3: J2000 RA (rad)
    # 4: J2000 DEC (rad) 5: horizontal az. (rad) 6: horizontal alt. (rad)
    # 7: source of test name
    testcases = [(48.86, 2.34, 1497003786, 0.8381070981, -0.5059006522, 3.26454625, 0.207694181, "stell. paris, alpha-for"),
                 (48.86, 2.34, 1497256975, 1.549728749, 0.1292784768, 1.969511946, 0.4898557422, "stell. paris, betelgeuse"),
                 (31.9633, -111.6, 2271646799, 0.05388704, 0.2650006125, 4.623697166, 0.66167856, "hor2eq IDL")]
    ARC_RAD = 0.00000484136811 # one arcsecond in rad
    for i in testcases:
        print 'Transformation test using:'
        print 'lat, long:', i[0], i[1], 'utc_time:', i[2], 'RA, DEC:', i[3], i[4]
        print 'az, alt:', i[5], i[6], 'from source: ', i[7]

        t = celestial.equatorial_to_horizontal_astropy(i[0], i[1], i[2], [(i[3], i[4])])
        print "equatorial_to_horizontal_astropy gave: ", t[0]
        print "Should have been:                      ", (i[5],i[6])
        print "Difference (dAZ, dALT) (arcsec):       ", ((t[0][0]-i[5])/ARC_RAD,(t[0][1]-i[6])/ARC_RAD)

        t = celestial.horizontal_to_equatorial_astropy(i[0], i[1], i[2], [(i[5], i[6])])
        print "horizontal_to_equatorial_astropy gave: ", t[0]
        print "Should have been:                      ", (i[3],i[4])
        print "Difference (dAZ, dALT) (arcsec):       ", ((t[0][0]-i[3])/ARC_RAD,(t[0][1]-i[4])/ARC_RAD)

        # equatorial to horizontal is actually a zenithazimuth to horizontal
        t = celestial.equatorial_to_horizontal(i[0], i[1], clock.utc_to_gps(i[2]), i[3], i[4])
        print "equatorial_to_horizontal gave:         ", (t[1], t[0])
        a,b = celestial.horizontal_to_zenithazimuth(i[6], i[5])
        print "Should have been:                      ", (b,a)
        print "Difference (dAZ, dALT) (arcsec):       ", ((t[1] - b) / ARC_RAD, (t[0] - a) / ARC_RAD)

        t = celestial.horizontal_to_equatorial(i[0], clock.utc_to_lst(datetime.datetime.utcfromtimestamp(i[2]), i[1]), i[6], i[5])
        print "horizontal_to_equatorial (is zenaz):   ", (t[0],t[1])
        print "Should have been:                      ", (i[3], i[4])
        print "Difference (dAZ, dALT) (arcsec):       ", ((t[0] - i[3]) / ARC_RAD, (t[1] - i[4]) / ARC_RAD)

        print("\n")

def oldvsnew_diagram():
    """
    Visual accuracy comparisons of old and new transformations.
    Compares the correlations between the transformations:
    equatorial_to_horizontal and equatorial_to_zenith_azimuth_astropy
    horizontal_to_equatorial and horizontal_to_zenith_azimuth_astropy
    Makes a histogram of the error differences between these two functions as well
    Assuming of course that astropy is relatively errorless and the old errorfull
    :return: None
    """
    # make random frames, in correct angle range and from utc time 2000-2020
    frames = []
    # boxes for the four different transformation results
    etoha = []
    etoh = []
    htoe = []
    htoea = []
    straight = lambda x : x # straight trendline function

    # Create the data sets for eq to az
    for i in range(100):
        frames.append((r.uniform(-90, 90), r.uniform(-180,180), r.randint(946684800,1577836800), r.uniform(0, 2*np.pi), r.uniform(-0.5*np.pi,0.5*np.pi)))
    for i in frames:
        etoha.append(celestial.equatorial_to_zenithazimuth_astropy(i[0],i[1], i[2], [(i[3], i[4])])[0])
        etoh.append(celestial.equatorial_to_horizontal(i[0], i[1], clock.utc_to_gps(i[2]), i[3], i[4]))
    # Data sets for hor to eq
    for i in frames:
        htoe.append(celestial.horizontal_to_equatorial(i[0], clock.utc_to_lst(datetime.datetime.utcfromtimestamp(i[2]), i[1]), i[4], i[3]))
        htoea.extend(celestial.horizontal_to_equatorial_astropy(i[0], i[1], i[2], [(i[3], i[4])]))

    # Make figs eq -> zenaz
    plt.figure(1)
    plt.suptitle('Zen/Az correlation in rads (equatorial_to_zenithazimuth/horizontal)')

    zenrange = [0, np.pi]
    plt.subplot(211)
    plt.title('Zenith')
    plt.axis(zenrange*2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')

    # Make figure and add 1:1 trendline

    plt.plot([co[0] for co in etoha], [co[0] for co in etoh], 'r.', zenrange, straight(zenrange), '-')

    plt.subplot(212)
    plt.title('Azimuth')
    azrange = [-np.pi, np.pi]
    plt.axis(azrange*2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')
    # Make figure and add 1:1 trendline
    plt.plot([co[1] for co in etoha], [co[1] for co in etoh], 'b.', azrange, straight(azrange), '-')
    plt.tight_layout() # Prevent titles merging
    plt.subplots_adjust(top=0.85)

    # Make histogram of differences
    plt.figure(2)
    nieuw = (np.array(etoh)-np.array(etoha))/2/np.pi*360*3600 # Take difference and convert to arcsec
    plt.hist([i[0] for i in nieuw], bins = 20)
    plt.title('Zenith Old-New Error (equatorial_to_zenithazimuth/horizontal)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    plt.figure(3)
    plt.hist([i[1] for i in nieuw], bins=20)
    plt.title('Azimuth Old-New Error (equatorial_to_zenithazimuth/horizontal)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    # Make figs hor - > eq

    plt.figure(4)
    plt.suptitle('RA/DEC correlation in rads (horizontal_to_equatorial)')
    altrange = [-0.5*np.pi, 0.5*np.pi]
    plt.subplot(211)
    plt.title('Declination')
    plt.axis(altrange * 2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')
    # Make figure and add 1:1 trendline
    plt.plot([co[1] for co in htoea], [co[1] for co in htoe], 'r.', altrange, straight(altrange), '-')

    plt.subplot(212)
    plt.title('Right Ascension')
    azrange = [0, 2*np.pi]
    plt.axis(azrange * 2)
    plt.xlabel('New (Astropy)')
    plt.ylabel('Old')
    # Make figure and add 1:1 trendline
    plt.plot([co[0] for co in htoea], [co[0] for co in htoe], 'b.', azrange, straight(azrange), '-')
    plt.tight_layout()  # Prevent titles merging
    plt.subplots_adjust(top = 0.85)

    # Make histogram of differences
    plt.figure(5)
    nieuw = (np.array(htoe) - np.array(htoea)) / 2 / np.pi * 360 * 3600  # Take difference and convert to arcsec
    plt.hist([i[1] for i in nieuw], bins=20)
    plt.title('Declination Old-New Error (horizontal_to_equatorial)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    plt.figure(6)
    nieuw = (np.array(htoe) - np.array(htoea)) / 2 / np.pi * 360 * 3600  # Take difference and convert to arcsec
    plt.hist([i[0] for i in nieuw], bins=20)
    plt.title('Right Ascension Old-New Error (horizontal_to_equatorial)')
    plt.xlabel('Error (arcsec)')
    plt.ylabel('Counts')

    plt.show()
    return

try:
    import ephem
    def pyephem_comp():
        # Set up randoms equatorial J2000 bodies that we will convert the RAs/Decs of.
        eq = [] # random frames to use
        for i in range(100):
            eq.append((r.uniform(-90, 90), r.uniform(-180, 180), r.randint(946684800, 1577836800),
                           r.uniform(0, 2 * np.pi), r.uniform(-0.5 * np.pi, 0.5 * np.pi)))
        efemeq = []
        altaz = []
        htoea = []
        etoha = []
        for i in eq:
            # Calculate altaz
            # Set observer for each case
            obs = ephem.Observer()
            obs.lat = str(i[0])
            obs.lon = str(i[1])
            obs.date = datetime.datetime.utcfromtimestamp(i[2])
            obs.pressure = 0 # Crucial to prevent refraction correction!

            # Set body for each case
            coord = ephem.FixedBody()
            coord._ra = i[3]
            coord._dec = i[4]

            # Do calculation and add to list, it is necessary to force the coords into radians
            coord.compute(obs)
            altaz.append((float(coord.az),float(coord.alt)))

            # Also calculate efemeq using eq
            result = obs.radec_of(i[3], i[4])
            efemeq.append((float(result[0]), float(result[1])))

        # Produce horizontal_to_equatorial_astropy results
        for i in eq:
            result = celestial.horizontal_to_equatorial_astropy(i[0], i[1], i[2], [(i[3],i[4])])
            htoea.extend(result)

        # Produce equatorial_to_horizontal_astropy results
        for i in eq:
            result = celestial.equatorial_to_horizontal_astropy(i[0], i[1], i[2], [(i[3], i[4])])
            etoha.extend(result)

        altdecrange = [-0.5*np.pi, 0.5*np.pi]
        azrarange = [0, 2*np.pi]

        plt.figure(1)
        plt.suptitle('RA/Dec correlation Altaz->(Astropy/Pyephem)->RA,DEC')

        # Create RA correlation subplot
        plt.subplot(211)
        plt.title('RA')
        plt.axis(azrarange*2)
        plt.xlabel('Pyephem RA (rad)')
        plt.ylabel('Astropy RA (rad)')

        # Plot RA correlation and trendline
        plt.plot([co[0] for co in efemeq], [co[0] for co in htoea], 'b.', azrarange, azrarange, '-')

        # DEC correlation subplot
        plt.subplot(212)
        plt.title('DEC')
        plt.axis(altdecrange*2)
        plt.xlabel('Pyephem DEC (rad)')
        plt.ylabel('Astropy DEC (rad)')

        # Plot DEC correlation and trendline
        plt.plot([co[1] for co in efemeq], [co[1] for co in htoea], 'r.', altdecrange, altdecrange, '-')

        # Formatting
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        # Make RA error histogram
        plt.figure(2)
        plt.title('RA Error Altaz->(Astropy/Pyephem)->RA,DEC')
        nieuw = (np.array(htoea) - np.array(efemeq))/2/np.pi* 360 * 3600 # Get differences in arcsec
        plt.hist([co[0] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')

        # Make DEC error histogram
        plt.figure(3)
        plt.title('DEC Error Altaz->(Astropy/Pyephem)->RA,DEC')
        plt.hist([co[1] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')

        # Altaz comparison plot
        plt.figure(4)
        plt.suptitle('Alt/Az correlation RA,DEC->(pyephem/astropy)->Altaz')

        # Altitude
        plt.subplot(211)
        plt.title('Altitude')
        plt.axis(altdecrange*2)
        plt.xlabel('Pyephem Altitude (rad)')
        plt.ylabel('Astropy Altitude (rad')

        # Plot with trendline
        plt.plot([co[1] for co in altaz], [co[1] for co in etoha], 'b.', altdecrange, altdecrange, '-')

        # Azimuth
        plt.subplot(212)
        plt.title('Azimuth')
        plt.axis(azrarange*2)
        plt.xlabel('Pyephem Azimuth (rad)')
        plt.ylabel('Astropy Azimuth (rad)')

        plt.plot([co[0] for co in altaz], [co[0] for co in etoha], 'r.', azrarange, azrarange, '-')

        # Formatting
        plt.tight_layout()
        plt.subplots_adjust(top=0.85)

        # Alt error histogram
        plt.figure(5)
        plt.title('Altitude Error RA,DEC->(pyephem/astropy)->Altaz')
        nieuw = (np.array(etoha)-np.array(altaz))/2/np.pi*360*3600
        plt.hist([co[1] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')

        # Az error histogram
        plt.figure(6)
        plt.title('Azimuth Error RA,DEC->(pyephem/astropy)->Altaz')
        plt.hist([co[0] for co in nieuw], bins=20)

        plt.ylabel('Counts')
        plt.xlabel('Error (arcsec)')
        # Done; output
        plt.show()


except ImportError:
    def pyephem_comp():
        print "Pyephem not present; no comparisons will be done"


if __name__ == '__main__':
    unittest.main()
