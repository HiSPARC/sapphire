import unittest
import random

import numpy as np

from sapphire.simulations import ldf


class BaseLdfSimulationTest(unittest.TestCase):

    def setUp(self):
        self.simulation = ldf.BaseLdfSimulation
        random.seed(1)
        np.random.seed(1)

    def test_simulate_particles_for_density(self):
        self.assertEqual(self.simulation.simulate_particles_for_density(100), 98)
        self.assertEqual(self.simulation.simulate_particles_for_density(4), 0)
        self.assertEqual(self.simulation.simulate_particles_for_density(4), 2)
        self.assertEqual(self.simulation.simulate_particles_for_density(0), 0)


class BaseLdfSimulationWithoutErrorsTest(BaseLdfSimulationTest):

    def setUp(self):
        super(BaseLdfSimulationWithoutErrorsTest, self).setUp()
        self.simulation = ldf.BaseLdfSimulationWithoutErrors

    def test_simulate_particles_for_density(self):
        self.assertEqual(self.simulation.simulate_particles_for_density(100), 100)
        self.assertEqual(self.simulation.simulate_particles_for_density(4), 4)
        self.assertEqual(self.simulation.simulate_particles_for_density(4), 4)
        self.assertEqual(self.simulation.simulate_particles_for_density(0), 0)


class BaseLdfTest(unittest.TestCase):

    def setUp(self):
        self.ldf = ldf.BaseLdf()

    def test_calculate_ldf_value(self):
        """Base LDF has no LDF, so no particles"""

        self.assertEqual(self.ldf.calculate_ldf_value(r=0), 0)
        self.assertEqual(self.ldf.calculate_ldf_value(r=10), 0)
        self.assertEqual(self.ldf.calculate_ldf_value(r=0, Ne=1e10), 0)
        self.assertEqual(self.ldf.calculate_ldf_value(r=0, Ne=1e10, s=3), 0)

    def test_calculate_core_distance(self):
        # TODO: Add core distances for inclined showers
        self.assertEqual(self.ldf.calculate_core_distance(0., 0., 0., 0., 0., 0.), 0)
        self.assertEqual(self.ldf.calculate_core_distance(10., 0., 0., 0., 0., 0.), 10.)
        self.assertEqual(self.ldf.calculate_core_distance(10., 0., 10., 0., 0., 0.), 0.)
        self.assertEqual(self.ldf.calculate_core_distance(10., 3., 10., 3., 0., 0.), 0.)


if __name__ == '__main__':
    unittest.main()
