import unittest
import tempfile
import os
import shutil
import tables
import sys
import cPickle as pickle

import numpy as np

from sapphire.analysis import direction_reconstruction


TEST_DATA_FILE = 'DIR-testdata.h5'
SIMULATION_GROUP = '/simulations/E_1PeV/zenith_22_5'


class DirectionReconstructionTests(unittest.TestCase):
    def setUp(self):
        self.data_path = self.create_tempfile_from_testdata()
        self.data = tables.open_file(self.data_path, 'a')

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
#        self.create_reconstruction_output(expected, overwrite=True)
        self.create_reconstruction_output(actual)

        self.validate_reconstruction_results(expected, actual)

    def test_binned_reconstruction_results(self):
        """Verify binned reconstruction results"""

        expected = '/reconstructions/prerecorded_binned'
        actual = '/reconstructions/test_binned'

        # For prerecording output, swap comments in following two lines
#        self.create_binned_reconstruction_output(expected, overwrite=True)
        self.create_binned_reconstruction_output(actual)

        self.validate_reconstruction_results(expected, actual)

    def test_randomized_binned_reconstruction_results(self):
        """Verify binned reconstruction results"""

        np.random.seed(1)

        expected = '/reconstructions/prerecorded_ran_binned'
        actual = '/reconstructions/test_ran_binned'

        # For prerecording output, swap comments in following two lines
#        self.create_randomized_binned_reconstruction_output(expected, overwrite=True)
        self.create_randomized_binned_reconstruction_output(actual)

        self.validate_reconstruction_results(expected, actual)

    def create_reconstruction_output(self, output, overwrite=False):
        reconstruction = direction_reconstruction.DirectionReconstruction(
            self.data, output, min_n134=1, N=100, overwrite=overwrite)
        self.redirect_stdout_stderr_to_devnull()
        reconstruction.reconstruct_angles_for_shower_group(SIMULATION_GROUP)
        self.restore_stdout_stderr()

    def create_binned_reconstruction_output(self, output, overwrite=False):
        reconstruction = direction_reconstruction.BinnedDirectionReconstruction(
            self.data, output, min_n134=1, N=100, binning=2.5, overwrite=overwrite)
        self.redirect_stdout_stderr_to_devnull()
        reconstruction.reconstruct_angles_for_shower_group(SIMULATION_GROUP)
        self.restore_stdout_stderr()

    def create_randomized_binned_reconstruction_output(self, output, overwrite=False):
        reconstruction = direction_reconstruction.BinnedDirectionReconstruction(
            self.data, output, min_n134=1, N=100, binning=2.5,
            randomize_binning=True, overwrite=overwrite)
        self.redirect_stdout_stderr_to_devnull()
        reconstruction.reconstruct_angles_for_shower_group(SIMULATION_GROUP)
        self.restore_stdout_stderr()

    def validate_reconstruction_results(self, expected, actual):
        expected = self.data.get_node(expected)
        actual = self.data.get_node(actual)
        self.validate_column_data(expected, actual)
        self.assertIn('cluster', actual.attrs)
        self.assertEqual(pickle.dumps(expected.attrs.cluster),
                         pickle.dumps(actual.attrs.cluster))

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
