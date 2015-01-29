import unittest
import tempfile
import os
import subprocess

import tables


TEST_DATA_PATH = 'test_data/'
TEST_EXPECTED_FILE = 'test_data/corsika_overview.h5'
STORE_SCRIPT = 'generate_corsika_overview {source} {destination}'


class GenerateCorsikaOverviewTests(unittest.TestCase):
    def setUp(self):
        self.source_path = self.get_testdata_path()
        self.expected_path = self.get_expected_path()
        self.destination_path = self.create_tempfile_path()
        self.command = STORE_SCRIPT.format(source=self.source_path,
                                           destination=self.destination_path)

    def tearDown(self):
        os.remove(self.destination_path)

    def test_store_data(self):
        result = subprocess.check_output(self.command, shell=True)
        self.assertEqual(result, '')
        self.validate_results()

    def validate_results(self):
        """Validate simulation results"""

        table_path = '/simulations'
        with tables.open_file(self.expected_path) as expected_file:
            with tables.open_file(self.destination_path) as actual_file:
                self.validate_table(table_path, expected_file, actual_file)

    def validate_table(self, table, expected_file, actual_file):
        """Verify that two tables are identical"""

        expected_node = expected_file.get_node(table)
        actual_node = actual_file.get_node(table)

        for colname in expected_node.colnames:
            expected_col = expected_node.col(colname)
            actual_col = actual_node.col(colname)
            if expected_col.shape == actual_col.shape:
                self.assertTrue((expected_col == actual_col).all())
            else:
                self.fail("Columns do not have the same length.")

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


if __name__ == '__main__':
    unittest.main()
