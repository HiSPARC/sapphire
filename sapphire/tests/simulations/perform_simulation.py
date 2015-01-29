"""Create test data for future acceptance testing"""

import os
import tempfile

import tables

import sapphire.clusters
from sapphire.simulations.groundparticles import GroundParticlesSimulation


self_path = os.path.dirname(__file__)
test_data_path = os.path.join(self_path, 'test_data/groundparticles_sim.h5')


def perform_simulation(filename):
    """Perform a small simulation and store results in filename"""

    corsika_data_path = os.path.join(self_path, 'test_data/corsika.h5')
    cluster = sapphire.clusters.SimpleCluster(size=40)
    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as datafile:
        sim = GroundParticlesSimulation(corsika_data_path, 70, cluster,
                                        datafile, N=10, seed=1)
        sim.run()


def create_tempfile_path():
    """Create a temporary file, close it, and return the path"""

    f, path = tempfile.mkstemp(suffix='.h5')
    os.close(f)
    return path


def create_and_store_test_data():
    """Create test data for future acceptance testing"""

    perform_simulation(test_data_path)


if __name__ == '__main__':
    create_and_store_test_data()
