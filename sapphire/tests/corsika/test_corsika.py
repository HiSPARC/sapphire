import unittest

from sapphire import corsika


DATA_FILE = 'DAT000000'


class CorsikaFileTests(unittest.TestCase):
    def setUp(self):
        self.file = corsika.CorsikaFile(DATA_FILE)

    def tearDown(self):
        pass

    def test_validate_file(self):
        """Verify that the data file is valid"""

        self.assertTrue(self.file.Check())

    def test_event(self):
        """Verify that the Event header is properly read"""

        events = self.file.GetEvents()
        event = events.next()
        header = event.GetHeader()
        self.assertEqual(header.fEnergy, 1e14)


if __name__ == '__main__':
    unittest.main()
