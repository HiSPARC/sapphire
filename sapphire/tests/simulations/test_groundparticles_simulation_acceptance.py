import os
import unittest
import warnings

from mock import patch

from sapphire.tests.validate_results import validate_results
from perform_simulation import (create_tempfile_path, perform_simulation,
                                test_data_path)


class GroundparticlesSimulationAcceptanceTest(unittest.TestCase):

    def setUp(self):
        warnings.filterwarnings('ignore')

    def tearDown(self):
        warnings.resetwarnings()

    @patch('sapphire.simulations.groundparticles.time')
    def test_simulation_output(self, mock_time):
        """Perform a simulation and verify the output"""

        mock_time.return_value = int(1e9)

        output_path = create_tempfile_path()
        perform_simulation(output_path)
        validate_results(self, test_data_path, output_path)
        os.remove(output_path)


if __name__ == '__main__':
    unittest.main()
