"""Create test data for future acceptance testing"""

import os
import random
import subprocess
import tempfile

import numpy
import tables

import sapphire.clusters
from sapphire.simulations.groundparticles import GroundParticlesSimulation


self_path = os.path.dirname(__file__)
test_data_path = os.path.join(self_path,
                              'test_data/groundparticles_sim.h5')


def perform_simulation(filename):
    """Perform a small simulation and store results in filename"""

    random.seed(1)
    numpy.random.seed(1)
    corsika_data_path = os.path.join(self_path, 'test_data/corsika.h5')
    cluster = sapphire.clusters.SimpleCluster(size=50)
    with tables.open_file(filename, 'w') as datafile:
        sim = GroundParticlesSimulation(corsika_data_path, 100, cluster,
                                        datafile, N=10)
        sim.run()


def create_tempfile_path():
    """Create a temporary file, close it, and return the path"""

    f, path = tempfile.mkstemp(suffix='.h5')
    os.close(f)
    return path


def create_and_store_test_data():
    """Create test data for future acceptance testing"""

    tmppath = create_tempfile_path()

    perform_simulation(tmppath)
    subprocess.check_call(['ptrepack', '-o', '--complevel', '1', tmppath,
                           test_data_path])


if __name__ == '__main__':
    create_and_store_test_data()
