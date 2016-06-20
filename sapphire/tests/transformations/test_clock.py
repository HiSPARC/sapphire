import unittest
import datetime
import random

from sapphire.transformations import clock


class DecimalTimeTests(unittest.TestCase):

    def test_time_to_decimal(self):
        """Check time to decimal hours conversion

        May be off by some microseconds due to precision

        """
        decimal = 11.51003771
        time_obj = datetime.time(11, 30, 36, 135756)

        self.assertEqual(clock.time_to_decimal(time_obj), decimal)
        self.assertEqual(clock.decimal_to_time(decimal), time_obj)
        self.assertEqual(clock.time_to_decimal(clock.decimal_to_time(decimal)), decimal)
        self.assertEqual(clock.decimal_to_time(clock.time_to_decimal(time_obj)), time_obj)


class DateTests(unittest.TestCase):

    def test_date_to_juliandate(self):
        self.assertEqual(clock.date_to_juliandate(2010, 12, 25), 2455555.5)

    def test_datetime_to_juliandate(self):
        self.assertEqual(clock.datetime_to_juliandate(datetime.datetime(2010, 12, 25, 12)),
                         2455556.0)


class ModifiedJulianDateTests(unittest.TestCase):

    def test_datetime_to_modifiedjd(self):
        mjd = clock.datetime_to_modifiedjd(datetime.datetime(2010, 12, 25, 12))
        self.assertEqual(mjd, 55555.5)

    def test_juliandate_to_modifiedjd(self):
        """Difference between Julian Date and Modified JD is 2400000.5"""

        self.assertEqual(clock.juliandate_to_modifiedjd(2400000.5), 0.)
        self.assertEqual(clock.modifiedjd_to_juliandate(0.), 2400000.5)

        for _ in range(5):
            modifiedjd = random.uniform(0, 5000000)
            self.assertAlmostEqual(
                clock.juliandate_to_modifiedjd(
                    clock.modifiedjd_to_juliandate(modifiedjd)),
                modifiedjd)
            juliandate = random.uniform(0, 5000000)
            self.assertAlmostEqual(
                clock.modifiedjd_to_juliandate(
                    clock.juliandate_to_modifiedjd(juliandate)),
                juliandate)


class JulianDateToDateTimeTests(unittest.TestCase):

    def test_juliandate_to_utc(self):
        self.assertEqual(clock.juliandate_to_utc(2400000.5),
                         datetime.datetime(1858, 11, 17))
        self.assertEqual(clock.juliandate_to_utc(2455581.40429),
                         datetime.datetime(2011, 1, 19, 21, 42, 10, 655997))

    def test_juliandate_to_utc_gap(self):
        self.assertEqual(clock.juliandate_to_utc(2299159.5),
                         datetime.datetime(1582, 10, 4))
        self.assertEqual(clock.juliandate_to_utc(2299160.5),
                         datetime.datetime(1582, 10, 15))

    def test_modifiedjd_to_utc(self):
        self.assertEqual(clock.modifiedjd_to_utc(55580.90429),
                         datetime.datetime(2011, 1, 19, 21, 42, 10, 655997))


class GMSTTests(unittest.TestCase):

    def test_utc_to_gmst(self):
        # Perhaps not perfect test, a few seconds of uncertainty exist..
        self.assertAlmostEqual(clock.utc_to_gmst(datetime.datetime(2010, 12, 25)),
                               clock.time_to_decimal(datetime.time(6, 13, 35, 852535)))


class LSTTests(unittest.TestCase):

    def test_gmst_to_lst(self):
        for _ in range(5):
            hours = random.uniform(0, 23.934)
            longitude = random.uniform(-180, 180)
            self.assertAlmostEqual(clock.lst_to_gmst(clock.gmst_to_lst(hours, longitude),
                                                     longitude), hours)

    def test_utc_to_lst_gmst(self):
        self.assertEqual(clock.utc_to_lst(datetime.datetime(2010, 12, 25), 0),
                         clock.utc_to_gmst(datetime.datetime(2010, 12, 25)))
        # Perhaps not perfect test, a few seconds of uncertainty exist..
        self.assertAlmostEqual(clock.utc_to_lst(datetime.datetime(2010, 12, 25), 0),
                               clock.time_to_decimal(datetime.time(6, 13, 35, 852535)))

    def test_utc_to_lst_at_longitudes(self):
        self.assertAlmostEqual(clock.utc_to_lst(datetime.datetime(2010, 12, 25), 90),
                               clock.time_to_decimal(datetime.time(12, 13, 35, 852535)))
        self.assertAlmostEqual(clock.utc_to_lst(datetime.datetime(2010, 12, 25), 180),
                               clock.time_to_decimal(datetime.time(18, 13, 35, 852535)))
        self.assertAlmostEqual(clock.utc_to_lst(datetime.datetime(2010, 12, 25), 5),
                               clock.time_to_decimal(datetime.time(6, 33, 35, 852535)))


class GPSTimeTests(unittest.TestCase):

    def setUp(self):
        """Setup combinations of calendar dates, timestamps and leap seconds

        The UTC timestamps follow from:

            import time
            import calendar
            t = calendar.timegm(time.strptime(date, '%B %d, %Y'))

        """
        self.combinations = (('July 1, 2015', 1435708800, 17),
                             ('January 1, 2014', 1388534400, 16),
                             ('July 1, 2012', 1341100800, 16),
                             ('June 30, 2012', 1341014400, 15),
                             ('January 1, 2009', 1230768000, 15),
                             ('December 31, 2008', 1230681600, 14),
                             ('January 1, 2006', 1136073600, 14),
                             ('December 31, 2005', 1135987200, 13),
                             ('January 1, 2004', 1072915200, 13),
                             ('January 1, 1999', 915148800, 13),
                             ('July 1, 1997', 867715200, 12),
                             ('January 1, 1996', 820454400, 11),
                             ('July 1, 1994', 773020800, 10),
                             ('July 1, 1993', 741484800, 9),
                             ('July 1, 1992', 709948800, 8),
                             ('January 1, 1991', 662688000, 7),
                             ('January 1, 1990', 631152000, 6),
                             ('January 1, 1988', 567993600, 5),
                             ('July 1, 1985', 489024000, 4),
                             ('July 1, 1983', 425865600, 3),
                             ('July 1, 1982', 394329600, 2),
                             ('July 1, 1981', 362793600, 1))

    def test_gps_to_utc(self):
        for date, _, _ in self.combinations:
            self.assertEqual(clock.gps_to_utc(clock.gps_from_string(date)),
                             clock.utc_from_string(date))

    def test_utc_to_gps(self):
        for date, _, _ in self.combinations:
            self.assertEqual(clock.utc_to_gps(clock.utc_from_string(date)),
                             clock.gps_from_string(date))

    def test_utc_from_string(self):
        for date, timestamp, _ in self.combinations:
            self.assertEqual(clock.utc_from_string(date), timestamp)

    def test_gps_from_string(self):
        for date, timestamp, leapseconds in self.combinations:
            self.assertEqual(clock.gps_from_string(date), timestamp + leapseconds)

    def test_gps_to_datetime(self):
        for date, timestamp, _ in self.combinations:
            dt = datetime.datetime.strptime(date, '%B %d, %Y')
            self.assertEqual(clock.gps_to_datetime(timestamp), dt)

    def test_datetime_to_gps(self):
        for date, timestamp, _ in self.combinations:
            dt = datetime.datetime.strptime(date, '%B %d, %Y')
            self.assertEqual(clock.datetime_to_gps(dt), timestamp)

    def test_process_time(self):
        for date, timestamp, _ in self.combinations:
            dt = datetime.datetime.strptime(date, '%B %d, %Y')
            self.assertEqual(clock.process_time(dt), timestamp)
            self.assertEqual(clock.process_time(timestamp), timestamp)
        self.assertEqual(clock.process_time('1435708800'), 1435708800)
        with self.assertRaises(RuntimeError):
            clock.process_time('July 1, 1995')


if __name__ == '__main__':
    unittest.main()
