from mock import sentinel, Mock, patch, MagicMock
import unittest

import tables

from sapphire.simulations.base import BaseSimulation
from sapphire import storage


class BaseSimulationTest(unittest.TestCase):
    def setUp(self):
        self.cluster = MagicMock(name='cluster')
        self.datafile = MagicMock(name='data')
        self.output_path = sentinel.output_path
        self.N = sentinel.N

        self.simulation = BaseSimulation(self.cluster, self.datafile,
                                         self.output_path, self.N)

    def test_init_sets_attributes(self):
        self.assertIs(self.simulation.cluster, self.cluster)
        self.assertIs(self.simulation.datafile, self.datafile)
        self.assertIs(self.simulation.output_path, self.output_path)
        self.assertIs(self.simulation.N, self.N)

    @patch.object(BaseSimulation, '_prepare_output_tables')
    def test_init_calls_prepare_output_tables(self, mock_method):
        BaseSimulation(Mock(), Mock(), Mock(), Mock())
        mock_method.assert_called_once_with()

    @unittest.skip("Does not test this unit")
    def test_init_creates_coincidences_output_group(self):
        self.datafile.createGroup.assert_any_call(self.output_path, 'coincidences', createparents=True)
        self.datafile.createTable.assert_called_with(self.simulation.coincidence_group, 'coincidences', storage.Coincidence)
        self.assertEqual(self.datafile.createVLArray.call_count, 2)
        self.datafile.createVLArray.assert_any_call(self.simulation.coincidence_group, 'c_index', tables.UInt32Col(shape=2))

    @unittest.skip("Does not test this unit")
    def test_init_creates_cluster_output_group(self):
        self.datafile.createGroup.assert_any_call(self.output_path, 'cluster_simulations', createparents=True)
        # The following tests need a better mock of cluster in order to work.
        # self.datafile.createGroup.assert_any_call(self.simulation.cluster_group, 'station_0')
        # self.datafile.createTable.assert_any_call(station_group, 'events', storage.ProcessedHisparcEvent, expectedrows=self.N)

    @unittest.skip("Does not test this unit")
    def test_init_stores_cluster_in_attrs(self):
        self.assertIs(self.simulation.coincidence_group._v_attrs.cluster, self.cluster)


if __name__ == '__main__':
    unittest.main()
