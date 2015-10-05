import unittest
import tempfile
import os
import subprocess

from sapphire.corsika.generate_corsika_overview import \
    generate_corsika_overview
from sapphire.tests.validate_results import validate_results

TEST_DATA_PATH = 'test_data/'
TEST_EXPECTED_FILE = 'test_data/corsika_overview.h5'
STORE_SCRIPT = 'generate_corsika_overview {source} {destination}'


class GenerateCorsikaOverviewTests(unittest.TestCase):

    def setUp(self):
        self.source_path = self.get_testdata_path()
        self.expected_path = self.get_expected_path()
        self.destination_path = self.create_tempfile_path()

    def tearDown(self):
        os.remove(self.destination_path)

    def test_store_data(self):
        generate_corsika_overview(source=self.source_path,
                                  destination=self.destination_path)
        validate_results(self, self.expected_path, self.destination_path)

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path

    def get_testdata_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_DATA_PATH)

    def get_expected_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_EXPECTED_FILE)


class GenerateCorsikaOverviewCommandTests(GenerateCorsikaOverviewTests):

    def setUp(self):
        self.source_path = self.get_testdata_path()
        self.expected_path = self.get_expected_path()
        self.destination_path = self.create_tempfile_path()
        self.command = STORE_SCRIPT.format(source=self.source_path,
                                           destination=self.destination_path)

    def test_store_data(self):
        result = subprocess.check_output(self.command, shell=True)
        self.assertEqual(result, '')
        validate_results(self, self.expected_path, self.destination_path)


if __name__ == '__main__':
    unittest.main()
