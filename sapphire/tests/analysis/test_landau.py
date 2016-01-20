import unittest

from numpy import linspace

from sapphire.analysis import landau


class LandauTest(unittest.TestCase):

    def test_pdf_mpv(self):
        """Check if peak of Landau pdf is at correct place

        The peak of the Landau should be around -0.22

        """
        x = linspace(-.4, 0, 100)
        self.assertAlmostEqual(x[landau.pdf(x).argmax()] + 0.222, 0, 2)


class ScintillatorTest(unittest.TestCase):

    def setUp(self):
        self.scin = landau.Scintillator()

    def test_pdf(self):
        """Check if the integral of the Landau pdf is almost 1"""

        self.scin.pdf(0)
        step_size = (self.scin.full_domain[-1] - self.scin.full_domain[-2])
        self.assertAlmostEqual(self.scin.pdf_values.sum() * step_size, 1, 1)


if __name__ == '__main__':
    unittest.main()
