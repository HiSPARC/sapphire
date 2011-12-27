import unittest
import tempfile
import os
import shutil
import tables

from numpy import deg2rad

import direction_reconstruction as rec


class DirectionReconstructionTests(unittest.TestCase):
    def test_reconstruction_results(self):
        """Verify that reconstruction output matches prerecorded output"""

        data_path = self.create_tempfile_from_testdata()
        self.data = tables.openFile(data_path, 'a')

        self.create_reconstruction_output()
        self.validate_reconstruction_results()

        self.data.close()
        os.remove(data_path)

    def create_reconstruction_output(self):
        output = self.data.createTable('/reconstructions',
                                       'test',
#                                       'prerecorded',
                                       rec.ReconstructedEvent, createparents=True)
        rec.reconstruct_angles(self.data, 'zenith_22_5', output, deg2rad(22.5), 1)

    def validate_reconstruction_results(self):
        expected = self.data.root.reconstructions.prerecorded
        actual = self.data.root.reconstructions.test
        self.validate_column_data(expected, actual)

    def validate_column_data(self, expected, actual):
        for colname in expected.colnames:
            expected_col = expected.col(colname)
            actual_col = actual.col(colname)
            self.assertTrue((expected_col == actual_col).all())

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
        return os.path.join(dir_path, 'DIR-testdata.h5')


if __name__ == '__main__':
    unittest.main()
