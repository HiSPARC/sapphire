from mock import Mock
import unittest
import os

import tables
from numpy import pi, sqrt

from sapphire.clusters import SingleDiamondStation
from sapphire.simulations.groundparticles import (GroundParticlesSimulation,
                                                  DetectorBoundarySimulation)


self_path = os.path.dirname(__file__)


class GroundParticlesSimulationTest(unittest.TestCase):

    def setUp(self):
        self.simulation = GroundParticlesSimulation.__new__(GroundParticlesSimulation)

        corsika_data_path = os.path.join(self_path, 'test_data/corsika.h5')
        self.corsika_data = tables.open_file(corsika_data_path, 'r')
        self.simulation.corsikafile = self.corsika_data

        self.simulation.cluster = SingleDiamondStation()
        self.detectors = self.simulation.cluster.stations[0].detectors

    def tearDown(self):
        self.corsika_data.close()

    def test__prepare_cluster_for_shower(self):

        # Combinations of shower parameters and detector after transformations
        combinations = (((0, 0, 0), (-0, -0, -0)),
                        ((10, -60, 0), (-10, 60, -0)),
                        ((10, -60, pi / 2), (60, 10, -pi / 2)))

        for input, expected in combinations:
            self.simulation._prepare_cluster_for_shower(*input)
            self.assertAlmostEqual(self.simulation.cluster.x, expected[0])
            self.assertAlmostEqual(self.simulation.cluster.y, expected[1])
            self.assertAlmostEqual(self.simulation.cluster.alpha, expected[2])

    def test_get_particles_query_string(self):
        self.simulation.groundparticles = Mock()

        # Combinations of shower parameters and detector after transformations
        combinations = (((0, 0, 0), (-0, -0, -0)),
                        ((10, -60, 0), (-10, 60, -0)),
                        ((10, -60, pi / 2), (60, 10, -pi / 2)))

        for input, expected in combinations:
            self.simulation._prepare_cluster_for_shower(*input)
            self.simulation.get_particles_in_detector(self.detectors[0])
            x, y = self.detectors[0].get_xy_coordinates()
            size = sqrt(.5) / 2.
            self.simulation.groundparticles.read_where.assert_called_with(
                '(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f) & '
                '(particle_id >= 2) & (particle_id <= 6)' %
                (x - size, x + size, y - size, y + size))

    def test_get_particles(self):
        self.groundparticles = self.corsika_data.root.groundparticles
        self.simulation.groundparticles = self.groundparticles

        combinations = (((0, 0, 0), (1, 0, 0, 0)),
                        ((1, -1, 0), (0, 1, 0, 3)),
                        ((1, -1, pi / 2), (1, 1, 0, 1)))

        for input, expected in combinations:
            self.simulation._prepare_cluster_for_shower(*input)
            for d, e in zip(self.detectors, expected):
                self.assertEqual(len(self.simulation.get_particles_in_detector(d)), e)


class DetectorBoundarySimulationTest(GroundParticlesSimulationTest):

    def setUp(self):
        self.simulation = DetectorBoundarySimulation.__new__(DetectorBoundarySimulation)

        corsika_data_path = os.path.join(self_path, 'test_data/corsika.h5')
        self.corsika_data = tables.open_file(corsika_data_path, 'r')
        self.simulation.corsikafile = self.corsika_data

        self.simulation.cluster = SingleDiamondStation()
        self.detectors = self.simulation.cluster.stations[0].detectors

    def test_get_particles_query_string(self):
        self.simulation.groundparticles = Mock()

        # Combinations of shower parameters and detector after transformations
        combinations = (((0, 0, 0), (-0, -0, -0)),
                        ((10, -60, 0), (-10, 60, -0)))

        for input, expected in combinations:
            self.simulation._prepare_cluster_for_shower(*input)
            self.simulation.get_particles_in_detector(self.detectors[0])
            x, y = self.detectors[0].get_xy_coordinates()
            size = .6
            self.simulation.groundparticles.read_where.assert_called_with(
                '(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f) & '
                '(b11 < y - 0.000000 * x) & (y - 0.000000 * x < b12) & '
                '(b21 < x) & (x < b22) & '
                '(particle_id >= 2) & (particle_id <= 6)' %
                (x - size, x + size, y - size, y + size))

    def test_get_particles(self):
        self.groundparticles = self.corsika_data.root.groundparticles
        self.simulation.groundparticles = self.groundparticles

        combinations = (((0, 0, 0), (1, 0, 1, 0)),
                        ((1, -1, 0), (0, 1, 1, 3)),
                        ((1, -1, pi / 2), (1, 1, 0, 1)))

        for input, expected in combinations:
            self.simulation._prepare_cluster_for_shower(*input)
            for d, e in zip(self.detectors, expected):
                self.assertEqual(len(self.simulation.get_particles_in_detector(d)), e)

    def test_get_line_boundary_eqs(self):
        combinations = ((((0, 0), (1, 1), (0, 2)), (0.0, 'y - 1.000000 * x', 2.0)),
                        (((0, 0), (0, 1), (1, 2)), (0.0, 'x', 1)))

        for input, expected in combinations:
            result = self.simulation.get_line_boundary_eqs(*input)
            self.assertEqual(result, expected)


if __name__ == '__main__':
    unittest.main()
