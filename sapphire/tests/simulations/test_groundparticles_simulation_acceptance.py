import os
import unittest

from mock import patch
import tables

from perform_simulation import (create_tempfile_path, perform_simulation,
                                test_data_path)


class GroundparticlesSimulationAcceptanceTest(unittest.TestCase):

    @patch('progressbar.ProgressBar')
    def test_simulation_output(self, mock_progressbar):
        """Perform a simulation and verify the output"""

        mock_progressbar.return_value.side_effect = lambda x: x

        output_path = create_tempfile_path()
        perform_simulation(output_path)
        self.validate_results(test_data_path, output_path)
        os.remove(output_path)

    def validate_results(self, expected_path, actual_path):
        """Validate simulation results"""

        expected_file = tables.open_file(expected_path)
        actual_file = tables.open_file(actual_path)
        
        for station_id in range(4):
            events_path = '/cluster_simulations/station_%d/events' % \
                          station_id
            self.validate_table(events_path, expected_file, actual_file)

        self.validate_table('/coincidences/coincidences', expected_file,
                            actual_file)
        self.validate_array('/coincidences/c_index', expected_file,
                            actual_file)
        self.validate_array('/coincidences/s_index', expected_file,
                            actual_file)

        expected_file.close()
        actual_file.close()

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

    def validate_array(self, table, expected_file, actual_file):
        """Verify that two arrays are identical"""

        expected_node = expected_file.get_node(table)
        actual_node = actual_file.get_node(table)

        expected = expected_node.read()
        actual = actual_node.read()

        self.assertEqual(len(expected), len(actual))

        if str(expected_node.atom) == 'VLStringAtom()':
            for i, j in zip(expected, actual):
                self.assertTrue(i == j)
        else:
            for i, j in zip(expected, actual):
                self.assertTrue((i == j).all())


if __name__ == '__main__':
    unittest.main()
