from mock import sentinel, Mock, patch, MagicMock
import unittest

from sapphire.old_simulations.base import BaseSimulation
import sapphire.clusters


class BaseSimulationTest(unittest.TestCase):
    @patch('os.path.split')
    def setUp(self, os_path_split_mock):
        self.output_head, self.output_tail = (sentinel.output_head,
                                              sentinel.output_tail)
        os_path_split_mock.return_value = self.output_head, self.output_tail
        self.output = sentinel.output

        self.cluster = sentinel.cluster
        self.data = MagicMock(name='data')
        self.R = sentinel.R
        self.N = sentinel.N

        self.simulation = BaseSimulation(self.cluster, self.data, self.output,
                                         self.R, self.N)

    def test_init_sets_attributes(self):
        self.assertIs(self.simulation.cluster, self.cluster)
        self.assertIs(self.simulation.data, self.data)
        self.assertIs(self.simulation.R, self.R)
        self.assertIs(self.simulation.N, self.N)

    def test_init_creates_output_group(self):
        self.data.createGroup.assert_called_with(self.output_head,
                                                 self.output_tail,
                                                 createparents=True)
        self.assertIs(self.simulation.output, self.data.createGroup.return_value)

    @patch('os.path.split')
    def test_init_raises_runtimeerror_if_output_exists(self, os_path_split_mock):
        os_path_split_mock.return_value = self.output_head, self.output_tail

        self.data.__contains__.return_value = True
        with self.assertRaises(RuntimeError):
            BaseSimulation(self.cluster, self.data, self.output, self.R, self.N)

    @patch('os.path.split')
    def test_init_remove_output_if_output_exists_and_force(self, os_path_split_mock):
        os_path_split_mock.return_value = self.output_head, self.output_tail

        self.data.__contains__.return_value = True
        BaseSimulation(self.cluster, self.data, self.output, self.R, self.N, force=True)
        self.data.removeNode.assert_called_with(self.output, recursive=True)

    def test_init_stores_cluster_in_attrs(self):
        self.assertIs(self.simulation.output._v_attrs.cluster, sentinel.cluster)


if __name__ == '__main__':
    unittest.main()
