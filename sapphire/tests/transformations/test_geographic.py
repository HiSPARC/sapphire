import unittest

from sapphire.transformations import geographic


class GeographicTransformationTests(unittest.TestCase):

    def setUp(self):
        self.ref_enu = (0., 0., 0.)
        self.ref_lla = (52.35592417, 4.95114402, 56.10234594)
        self.transform = geographic.FromWGS84ToENUTransformation(self.ref_lla)

    def test_attributes(self):
        self.assertEqual(self.ref_lla, self.transform.ref_lla)

        ref_ecef = (3889144.77, 336914.68, 5027133.30)
        self.assert_tuple_almost_equal(ref_ecef, self.transform.ref_ecef, 2)

    def test_lla_to_enu_at_reference(self):
        self.assert_tuple_almost_equal(self.ref_enu, self.transform.lla_to_enu(self.transform.ref_lla))

    def test_enu_to_lla_at_reference(self):
        self.assert_tuple_almost_equal(self.ref_lla, self.transform.enu_to_lla(self.ref_enu))

    def test_lla_to_ecef_at_reference(self):
        self.assert_tuple_almost_equal(self.transform.ref_ecef, self.transform.lla_to_ecef(self.ref_lla))

    def test_ecef_to_lla_at_reference(self):
        self.assert_tuple_almost_equal(self.ref_lla, self.transform.ecef_to_lla(self.transform.ref_ecef))

    def test_ecef_to_enu_at_reference(self):
        self.assert_tuple_almost_equal(self.ref_enu, self.transform.ecef_to_enu(self.transform.ref_ecef))

    def test_enu_to_ecef_at_reference(self):
        self.assert_tuple_almost_equal(self.transform.ref_ecef, self.transform.enu_to_ecef(self.ref_enu))

    def assert_tuple_almost_equal(self, actual, expected, places=7):
        self.assertIsInstance(actual, tuple)
        self.assertIsInstance(expected, tuple)

        msg = f"Tuples differ: {actual} != {expected}"
        for actual_value, expected_value in zip(actual, expected):
            self.assertAlmostEqual(actual_value, expected_value, places=places, msg=msg)


if __name__ == '__main__':
    unittest.main()
