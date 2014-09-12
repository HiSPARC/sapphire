import unittest

from numpy import pi

from sapphire.transformations import angles


class DecimalSexagesimalTests(unittest.TestCase):

    def setUp(self):
        self.combinations = ((0, (0, 0, 0)),
                             (1, (1, 0, 0)),
                             (30, (30, 0, 0)),
                             (1 / 60., (0, 1, 0)),
                             (.5, (0, 30, 0)),
                             (1 / 3600., (0, 0, 1)),
                             (30 / 3600., (0, 0, 30)))

    def test_decimal_to_sexagesimal(self):
        for dec, sexa in self.combinations:
            self.assertEqual(angles.decimal_to_sexagesimal(dec), sexa)

    def test_sexagesimal_to_decimal(self):
        for dec, sexa in self.combinations:
            self.assertEqual(angles.sexagesimal_to_decimal(*sexa), dec)


if __name__ == '__main__':
    unittest.main()
