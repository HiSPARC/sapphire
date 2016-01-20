import unittest
import tempfile
import os

import tables

from sapphire import kascade
from sapphire.tests.validate_results import validate_results


self_path = os.path.dirname(__file__)
TEST_DATA_FILE = os.path.join(self_path, 'test_data/kascade.dat')
TEST_DATA_REF = os.path.join(self_path, 'test_data/kascade.h5')


class StoreKascadeDataTests(unittest.TestCase):

    def setUp(self):
        self.destination_path = self.create_tempfile_path()

    def test_read_and_store_data(self):
        path = self.destination_path
        with tables.open_file(path, 'a') as self.destination_data:
            self.kascade = kascade.StoreKascadeData(self.destination_data,
                                                    TEST_DATA_FILE, '/kascade',
                                                    progress=False)
            self.kascade.read_and_store_data()
        validate_results(self, TEST_DATA_REF, self.destination_path)

    def tearDown(self):
        os.remove(self.destination_path)

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path


if __name__ == '__main__':
    unittest.main()
