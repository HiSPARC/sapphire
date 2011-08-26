import unittest

import tables
from mock import Mock, patch

from simulations import base

class BaseSimulationTests(unittest.TestCase):
    @patch('os.path.split')
    def setUp(self, os_path_split_mock):
        self.output_head, self.output_tail = (Mock(name='output_head'),
                                              Mock(name='output_tail'))
        os_path_split_mock.return_value = self.output_head, self.output_tail

        self.cluster = Mock(name='cluster')
        self.data = Mock(name='data')
        self.grdpcles = Mock(name='grdpcles')
        self.output = Mock(name='output')
        self.R = Mock(name='R')
        self.N = Mock(name='N')
        self.simulation = base.BaseSimulation(self.cluster, self.data,
                                              self.grdpcles, self.output,
                                              self.R, self.N)

    def test_init_sets_attributes(self):
        self.assertIs(self.simulation.cluster, self.cluster)
        self.assertIs(self.simulation.data, self.data)
        self.assertIs(self.simulation.R, self.R)
        self.assertIs(self.simulation.N, self.N)

    def test_init_gets_grdpcles_node(self):
        self.data.getNode.assert_called_with('/', self.grdpcles)
        self.assertIs(self.simulation.grdpcles, self.data.getNode.return_value)

    @patch('os.path.split')
    def test_init_raises_runtimeerror_if_grdpcles_not_found(self, os_path_split_mock):
        os_path_split_mock.return_value = Mock(), Mock()

        data = Mock()
        data.getNode.side_effect = tables.NoSuchNodeError
        with self.assertRaises(RuntimeError):
            base.BaseSimulation(Mock(), data, Mock(), Mock(), Mock(), Mock())

    def test_init_creates_output(self):
        self.data.createGroup.assert_called_with(self.output_head,
                                                 self.output_tail,
                                                 createparents=True)
        self.assertIs(self.simulation.output, self.data.createGroup.return_value)

    @patch('os.path.split')
    def test_init_raises_runtimeerror_if_creategroup_error(self, os_path_split_mock):
        os_path_split_mock.return_value = Mock(), Mock()

        data = Mock()
        data.createGroup.side_effect = tables.NodeError
        with self.assertRaises(RuntimeError):
            base.BaseSimulation(Mock(), data, Mock(), Mock(), Mock(), Mock())
