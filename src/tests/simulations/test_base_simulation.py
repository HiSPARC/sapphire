from mock import sentinel, Mock, patch, MagicMock
import unittest

from simulations.base import BaseSimulation
import clusters


class BaseSimulationTest(unittest.TestCase):
    def setUp(self):
        cluster = sentinel.cluster
        data = MagicMock()
        output = '/simulations'
        R = sentinel.R
        N = sentinel.N

        self.simulation = BaseSimulation(cluster, data, output, R, N)

    def test_init_stores_cluster_in_attrs(self):
        self.assertIs(self.simulation.output._v_attrs.cluster, sentinel.cluster)
