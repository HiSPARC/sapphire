import unittest
import tempfile
import os
import shutil
import tables
import sys

from numpy import deg2rad
import numpy as np

import direction_reconstruction


TEST_DATA_FILE = 'DIR-testdata.h5'
SIMULATION_GROUP = '/simulations/E_1PeV/zenith_22_5'


class DirectionReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.openFile(self.data_path, 'a')

    def tearDown(self):
        self.data.close()

        # For prerecording output, swap comments in following two lines
#        os.rename(self.data_path, self.get_testdata_path())
        os.remove(self.data_path)

    def test_reconstruction_results(self):
        """Verify that reconstruction output matches prerecorded output"""

        expected = '/reconstructions/prerecorded'
        actual = '/reconstructions/test'

        # For prerecording output, swap comments in following two lines
#        self.create_reconstruction_output(expected)
        self.create_reconstruction_output(actual)

        self.validate_reconstruction_results(expected, actual)

    def test_binned_reconstruction_results(self):
        """Verify binned reconstruction results"""

        expected = '/reconstructions/prerecorded_binned'
        actual = '/reconstructions/test_binned'

        # For prerecording output, swap comments in following two lines
#        self.create_binned_reconstruction_output(expected)
        self.create_binned_reconstruction_output(actual)

        self.validate_reconstruction_results(expected, actual)

    def test_randomized_binned_reconstruction_results(self):
        """Verify binned reconstruction results"""

        np.random.seed(1)

        expected = '/reconstructions/prerecorded_ran_binned'
        actual = '/reconstructions/test_ran_binned'

        # For prerecording output, swap comments in following two lines
#        self.create_randomized_binned_reconstruction_output(expected)
        self.create_randomized_binned_reconstruction_output(actual)

        self.validate_reconstruction_results(expected, actual)

    def create_reconstruction_output(self, table_path):
        output = self.create_empty_output_table(table_path)
        self.reconstruct_direction(output)

    def create_binned_reconstruction_output(self, table_path):
        output = self.create_empty_output_table(table_path)
        self.reconstruct_direction_binned(output)

    def create_randomized_binned_reconstruction_output(self, table_path):
        output = self.create_empty_output_table(table_path)
        self.reconstruct_direction_randomized_binned(output)

    def create_empty_output_table(self, table_path):
        group, tablename = os.path.split(table_path)

        if not group in self.data.root:
            base, groupname = os.path.split(group)
            self.data.createGroup(base, groupname, createparents=True)
        group = self.data.getNode(group)

        if tablename in group:
            self.data.removeNode(table_path)

        output = self.data.createTable(group, tablename, direction_reconstruction.ReconstructedEvent)
        return output

    def reconstruct_direction(self, output):
        reconstruction = direction_reconstruction.DirectionReconstruction(
                            self.data, output, min_n134=1, N=100)
        self.redirect_stdout_stderr_to_devnull()
        reconstruction.reconstruct_angles_for_group(SIMULATION_GROUP, deg2rad(22.5))
        self.restore_stdout_stderr()

    def reconstruct_direction_binned(self, output):
        reconstruction = direction_reconstruction.DirectionReconstruction(
                            self.data, output, min_n134=1, N=100)
        self.redirect_stdout_stderr_to_devnull()
        reconstruction.reconstruct_angles_for_group(SIMULATION_GROUP, deg2rad(22.5), binning=2.5)
        self.restore_stdout_stderr()

    def reconstruct_direction_randomized_binned(self, output):
        reconstruction = direction_reconstruction.DirectionReconstruction(
                            self.data, output, min_n134=1, N=100)
        self.redirect_stdout_stderr_to_devnull()
        reconstruction.reconstruct_angles_for_group(SIMULATION_GROUP, deg2rad(22.5), binning=2.5, randomize_binning=True)
        self.restore_stdout_stderr()

    def validate_reconstruction_results(self, expected, actual):
        expected = self.data.getNode(expected)
        actual = self.data.getNode(actual)
        self.validate_column_data(expected, actual)

    def validate_column_data(self, expected, actual):
        for colname in expected.colnames:
            expected_col = expected.col(colname)
            actual_col = actual.col(colname)
            self.assertTrue((expected_col == actual_col).all())
            self.assertIs(expected_col.dtype, actual_col.dtype)

    def create_tempfile_from_testdata(self):
        tmp_path = self.create_tempfile_path()
        data_path = self.get_testdata_path()
        shutil.copyfile(data_path, tmp_path)
        return tmp_path

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path

    def get_testdata_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_DATA_FILE)

    def redirect_stdout_stderr_to_devnull(self):
        self.__stdout = sys.stdout
        self.__stderr = sys.stderr
        sys.stdout = open(os.devnull, 'w')
        sys.stderr = open(os.devnull, 'w')

    def restore_stdout_stderr(self):
        sys.stdout.close()
        sys.stderr.close()
        sys.stdout = self.__stdout
        sys.stderr = self.__stderr


if __name__ == '__main__':
    unittest.main()
