import unittest
import tables
import tempfile
import os
import numpy as np
import sys

from sapphire.old_simulations.ldf import BaseLdfSimulation
from sapphire import clusters


class LdfSimulationAcceptanceTest(unittest.TestCase):
    def test_simulation_yields_results(self):
        data_path = self.create_tempfile_path()
        self.data = tables.open_file(data_path, 'w')

        self.create_test_simulation_output()

        self.data.close()
        os.remove(data_path)

    def create_tempfile_path(self):
        fd, path = tempfile.mkstemp('.h5')
        os.close(fd)
        return path

    def create_test_simulation_output(self):
        np.random.seed(1)
        cluster = clusters.SimpleCluster()
        simulation = BaseLdfSimulation(cluster, self.data, '/sim', R=100, N=100)
        self.redirect_stdout_stderr_to_devnull()
        simulation.run()
        self.restore_stdout_stderr()

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
