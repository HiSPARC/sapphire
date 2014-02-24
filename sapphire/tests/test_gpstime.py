import unittest

from sapphire import gpstime


class GPSTimeTests(unittest.TestCase):
    def setUp(self):
        """Setup combinations of calendar dates, timestamps and leap seconds

        The UTC timestamps follow from:

            import time
            import calendar
            t = calendar.timegm(time.strptime(date, '%B %d, %Y'))

        """
        self.combinations = (('January 1, 1999', 915148800, 13),
                             ('January 1, 2004', 1072915200, 13),
                             ('December 31, 2005', 1135987200, 13),
                             ('January 1, 2006', 1136073600, 14),
                             ('December 31, 2008', 1230681600, 14),
                             ('January 1, 2009', 1230768000, 15),
                             ('June 30, 2012', 1341014400, 15),
                             ('July 1, 2012', 1341100800, 16),
                             ('January 1, 2014', 1388534400, 16))

    def test_out_of_range(self):
        self.assertRaises(Exception, gpstime.gps_to_utc,
                          gpstime.gps_from_string('January 1, 1999') - 1)
        self.assertRaises(Exception, gpstime.utc_to_gps,
                          gpstime.utc_from_string('January 1, 1999') - 1)

    def test_gps_to_utc(self):
        for date, _, _ in self.combinations:
            self.assertEqual(gpstime.gps_to_utc(gpstime.gps_from_string(date)),
                             gpstime.utc_from_string(date))

    def test_utc_to_gps(self):
        for date, _, _ in self.combinations:
            self.assertEqual(gpstime.utc_to_gps(gpstime.utc_from_string(date)),
                             gpstime.gps_from_string(date))

    def test_utc_from_string(self):
        for date, timestamp, _ in self.combinations:
            self.assertEqual(gpstime.utc_from_string(date), timestamp)

    def test_gps_from_string(self):
        for date, timestamp, leapseconds in self.combinations:
            self.assertEqual(gpstime.gps_from_string(date), timestamp + leapseconds)


if __name__ == '__main__':
    unittest.main()
