import unittest
import tempfile
import os
import subprocess

from sapphire.corsika.store_corsika_data import store_and_sort_corsika_data
from sapphire.tests.validate_results import validate_results

TEST_DATA_FILE = 'test_data/1_2/DAT000000'
TEST_EXPECTED_FILE = 'test_data/1_2/corsika.h5'
STORE_CMD = 'store_corsika_data {source} {destination}'


class StoreCorsikaDataTests(unittest.TestCase):

    """Store CORSIKA test using the function directly"""

    def setUp(self):
        self.source_path = self.get_testdata_path()
        self.expected_path = self.get_expected_path()
        self.destination_path = self.create_tempfile_path()

    def tearDown(self):
        os.remove(self.destination_path)

    def test_store_data(self):
        # First with overwrite false
        self.assertRaises(Exception, store_and_sort_corsika_data,
                          self.source_path, self.destination_path,
                          progress=True)
        # Now with overwrite true
        store_and_sort_corsika_data(self.source_path, self.destination_path,
                                    overwrite=True)
        validate_results(self, self.expected_path, self.destination_path)

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path

    def get_testdata_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_DATA_FILE)

    def get_expected_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_EXPECTED_FILE)


class StoreCorsikaDataCommandTests(StoreCorsikaDataTests):

    """Store CORSIKA test calling store command"""

    def setUp(self):
        self.source_path = self.get_testdata_path()
        self.expected_path = self.get_expected_path()
        self.destination_path = self.create_tempfile_path()
        self.command = STORE_CMD.format(source=self.source_path,
                                        destination=self.destination_path)

    def test_store_data(self):
        result = subprocess.check_output(self.command, shell=True)
        self.assertEqual(result, '')

        self.assertRaises(subprocess.CalledProcessError,
                          subprocess.check_output,
                          self.command + ' --progress',
                          stderr=subprocess.STDOUT, shell=True)

        result = subprocess.check_output(self.command + ' --overwrite',
                                         shell=True)
        self.assertEqual(result, '')

        validate_results(self, self.expected_path, self.destination_path)


if __name__ == '__main__':
    unittest.main()
