import unittest
import datetime
import random

from sapphire.transformations import time


class DecimalTimeTests(unittest.TestCase):

    def test_time_to_decimal(self):
        """Check time to decimal hours conversion

        May be off by some microseconds due to precision

        """
        decimal = 11.51003771
        time_obj = datetime.time(11, 30, 36, 135756)

        self.assertEqual(time.time_to_decimal(time_obj), decimal)
        self.assertEqual(time.decimal_to_time(decimal), time_obj)
        self.assertEqual(time.time_to_decimal(time.decimal_to_time(decimal)), decimal)
        self.assertEqual(time.decimal_to_time(time.time_to_decimal(time_obj)), time_obj)


class DateTests(unittest.TestCase):

    def test_date_to_juliandate(self):
        self.assertEqual(time.date_to_juliandate(2010, 12, 25), 2455555.5)

    def test_datetime_to_juliandate(self):
        self.assertEqual(time.datetime_to_juliandate(datetime.datetime(2010, 12, 25, 12)),
                         2455556.0)

class ModifiedJulianDateTests(unittest.TestCase):

    def test_juliandate_to_modifiedjd(self):
        """Difference between Julian Date and Modified JD is 2400000.5"""

        self.assertEqual(time.juliandate_to_modifiedjd(2400000.5), 0.)
        self.assertEqual(time.modifiedjd_to_juliandate(0.), 2400000.5)

        modifiedjd = random.uniform(0, 5000000)
        self.assertAlmostEqual(time.juliandate_to_modifiedjd(
                                   time.modifiedjd_to_juliandate(modifiedjd)),
                               modifiedjd)
        juliandate = random.uniform(0, 5000000)
        self.assertAlmostEqual(time.modifiedjd_to_juliandate(
                                   time.juliandate_to_modifiedjd(juliandate)),
                               juliandate)


class JulianDateToDateTimeTests(unittest.TestCase):

    def test_juliandate_to_utc(self):
        self.assertEqual(time.juliandate_to_utc(2400000.5),
                         datetime.datetime(1858, 11, 17))
        self.assertEqual(time.juliandate_to_utc(2455581.40429),
                         datetime.datetime(2011, 1, 19, 21, 42, 10, 655997))

    def test_modifiedjd_to_utc(self):
        self.assertEqual(time.modifiedjd_to_utc(55580.90429),
                         datetime.datetime(2011, 1, 19, 21, 42, 10, 655997))


class GMSTTests(unittest.TestCase):

    def test_utc_to_gmst(self):
        # Perhaps not perfect test, a few seconds of uncertainty exist..
        self.assertAlmostEqual(time.utc_to_gmst(datetime.datetime(2010, 12, 25)),
                               time.time_to_decimal(datetime.time(6, 13, 35, 852535)))

class LSTTests(unittest.TestCase):

    def test_utc_to_lst_gmst(self):
        self.assertEqual(time.utc_to_lst(0, datetime.datetime(2010, 12, 25)),
                         time.utc_to_gmst(datetime.datetime(2010, 12, 25)))
        # Perhaps not perfect test, a few seconds of uncertainty exist..
        self.assertAlmostEqual(time.utc_to_lst(0, datetime.datetime(2010, 12, 25)),
                               time.time_to_decimal(datetime.time(6, 13, 35, 852535)))

    def test_utc_to_lst_at_longitudes(self):
        self.assertAlmostEqual(time.utc_to_lst(90, datetime.datetime(2010, 12, 25)),
                               time.time_to_decimal(datetime.time(12, 13, 35, 852535)))
        self.assertAlmostEqual(time.utc_to_lst(180, datetime.datetime(2010, 12, 25)),
                               time.time_to_decimal(datetime.time(18, 13, 35, 852535)))
        self.assertAlmostEqual(time.utc_to_lst(5, datetime.datetime(2010, 12, 25)),
                               time.time_to_decimal(datetime.time(6, 33, 35, 852535)))



if __name__ == '__main__':
    unittest.main()
