import unittest
import tempfile
import os
import subprocess

from mock import patch


TEST_DATA_FILE = 'DAT000000'
STORE_SCRIPT = 'store_corsika_data {source} {destination}'

@unittest.skipIf(not os.path.exists(os.path.join(os.path.dirname(__file__), TEST_DATA_FILE)),
                 'Missing test datafile.')
class StoreCorsikaDataTests(unittest.TestCase):
    def setUp(self):
        self.source_path = self.get_testdata_path()
        self.destination_path = self.create_tempfile_path()
        self.command = STORE_SCRIPT.format(source=self.source_path,
                                           destination=self.destination_path)

    def tearDown(self):
        os.remove(self.destination_path)

    def test_store_data(self):
        result = subprocess.check_output(self.command, shell=True)
        self.assertEqual(result, '')

        self.assertRaises(subprocess.CalledProcessError,
                          subprocess.check_output, self.command,
                          stderr=subprocess.STDOUT, shell=True)

        result = subprocess.check_output(self.command + ' --overwrite', shell=True)
        self.assertEqual(result, '')

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path

    def get_testdata_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_DATA_FILE)


if __name__ == '__main__':
    unittest.main()
