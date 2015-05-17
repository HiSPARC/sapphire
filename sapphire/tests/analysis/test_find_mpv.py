import unittest
import warnings

from numpy import array

from sapphire.analysis import find_mpv


class FindMostProbableValueInSpectrumTest(unittest.TestCase):

    def test_failing_fit(self):
        """Check for correct warnings/errors for failing fit"""

        n = array([0, 1, 1])
        bins = array([0, 10, 20])
        fmpv = find_mpv.FindMostProbableValueInSpectrum(n, bins)

        # Exception from the fit mpv function, to few points
        first_guess = fmpv.find_first_guess_mpv()
        with self.assertRaises(RuntimeError) as cm:
            fmpv.fit_mpv(first_guess)
        self.assertEqual(str(cm.exception),
                         "Number of data points not sufficient")

        # Warning from the find mpv function
        with warnings.catch_warnings(record=True) as w:
            # clear the warnings from find_mpv module
            # http://bugs.python.org/issue4180
            if hasattr(find_mpv, '__warningregistry__'):
                find_mpv.__warningregistry__ = {}
            warnings.simplefilter("always")
            mpv, is_fitted = fmpv.find_mpv()
        self.assertTrue(issubclass(w[0].category, UserWarning))
        self.assertEqual(mpv, -999)
        self.assertFalse(is_fitted)

    @unittest.skip('Need better test, this has different error on Travis.')
    def test_bad_fit(self):
        """Exception from the fit mpv function, result outside range"""

        n = array([1, 3, 7, 70])
        bins = array([111.0, 111.1, 111.2, 111.3])
        fmpv = find_mpv.FindMostProbableValueInSpectrum(n, bins)
        with self.assertRaises(RuntimeError) as cm:
            fmpv.fit_mpv(111)
        self.assertEqual(str(cm.exception), "Fitted MPV value outside range")


if __name__ == '__main__':
    unittest.main()
