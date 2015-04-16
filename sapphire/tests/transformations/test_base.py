import unittest

from sapphire.transformations import base


class DecimalSexagesimalTests(unittest.TestCase):

    def setUp(self):
        # (decimal, sexagesimal)
        self.combinations = ((0, (0, 0, 0)),
                             (1, (1, 0, 0)),
                             (30, (30, 0, 0)),
                             (1 / 60., (0, 1, 0)),
                             (-1 + (30 / 60.), (0, -30, 0)),
                             (-1 - (30 / 60.) - (30 / 3600.), (-1, -30, -30)),
                             (.5, (0, 30, 0)),
                             (1 / 3600., (0, 0, 1)),
                             (30 / 3600., (0, 0, 30)))

    def test_decimal_to_sexagesimal(self):
        for dec, sexa in self.combinations:
            self.assertEqual(base.decimal_to_sexagesimal(dec), sexa)

    def test_sexagesimal_to_decimal(self):
        for dec, sexa in self.combinations:
            self.assertEqual(base.sexagesimal_to_decimal(*sexa), dec)


if __name__ == '__main__':
    unittest.main()
