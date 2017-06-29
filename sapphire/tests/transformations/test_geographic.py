import unittest

from sapphire.transformations import geographic


class GeographicTransformationTests(unittest.TestCase):

    def setUp(self):
        self.ref_lla = (52.35592417, 4.95114402, 56.10234594)
        self.transform = geographic.FromWGS84ToENUTransformation(self.ref_lla)

    def test_attributes(self):
        self.assertEqual(self.ref_lla, self.transform.ref_lla)

        ref_ecef = (3889144.77, 336914.68, 5027133.30)
        self.assert_tuple_almost_equal(ref_ecef, self.transform.ref_ecef, 2)

    def test_ecef_to_lla(self):
        lla = self.transform.ecef_to_lla(self.transform.ref_ecef)
        self.assert_tuple_almost_equal(self.ref_lla, lla)

    def test_ecef_to_enu(self):
        enu = self.transform.ecef_to_enu(self.transform.ref_ecef)
        ecef = self.transform.enu_to_ecef(enu)
        self.assert_tuple_almost_equal(self.transform.ref_ecef, ecef)

    def test_lla_to_enu(self):
        enu = self.transform.lla_to_enu(self.transform.ref_lla)
        lla = self.transform.enu_to_lla(enu)
        self.assert_tuple_almost_equal(self.transform.ref_lla, lla)

    def assert_tuple_almost_equal(self, actual, expected, places=7):
        self.assertIsInstance(actual, tuple)
        self.assertIsInstance(expected, tuple)

        msg = "Tuples differ: %s != %s" % (str(actual), str(expected))
        for actual_value, expected_value in zip(actual, expected):
            self.assertAlmostEqual(actual_value, expected_value, places=places,
                                   msg=msg)


if __name__ == '__main__':
    unittest.main()
