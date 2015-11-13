"""Create test data for future acceptance testing"""

import os
import tempfile

import tables
from mock import patch

import sapphire.clusters
from sapphire.simulations.groundparticles import GroundParticlesSimulation
from sapphire.simulations.showerfront import FlatFrontSimulation
from sapphire.simulations.ldf import NkgLdfSimulation


self_path = os.path.dirname(__file__)
test_data_path = os.path.join(self_path, 'test_data/groundparticles_sim.h5')
test_data_flat = os.path.join(self_path, 'test_data/flatfront_sim.h5')
test_data_nkg = os.path.join(self_path, 'test_data/nkgldf_sim.h5')


@patch('sapphire.simulations.groundparticles.time')
def perform_groundparticlessimulation(filename, mock_time):
    """Perform a small simulation and store results in filename"""

    mock_time.return_value = int(1e9)

    corsika_data_path = os.path.join(self_path, 'test_data/corsika.h5')
    cluster = sapphire.clusters.SimpleCluster(size=40)
    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as data:
        sim = GroundParticlesSimulation(corsika_data_path, 70, cluster,
                                        data, N=10, seed=1, progress=False)
        sim.run()


def perform_flatfrontsimulation(filename):
    """Perform a small simulation and store results in filename"""

    cluster = sapphire.clusters.SimpleCluster(size=40)
    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as data:
        sim = FlatFrontSimulation(cluster, data, '/', 10, seed=1,
                                  progress=False)
        sim.run()


def perform_nkgldfsimulation(filename):
    """Perform a small simulation and store results in filename"""

    cluster = sapphire.clusters.SimpleCluster(size=40)
    filters = tables.Filters(complevel=1)
    with tables.open_file(filename, 'w', filters=filters) as data:
        sim = NkgLdfSimulation(400, 1e15, 1e19, cluster, data, '/', 10,
                               seed=1, progress=False)
        sim.run()


def create_tempfile_path():
    """Create a temporary file, close it, and return the path"""

    f, path = tempfile.mkstemp(suffix='.h5')
    os.close(f)
    return path


def create_and_store_test_data():
    """Create reference test data for future acceptance testing"""

    perform_groundparticlessimulation(test_data_path)
    perform_flatfrontsimulation(test_data_flat)
    perform_nkgldfsimulation(test_data_nkg)


if __name__ == '__main__':
    create_and_store_test_data()
