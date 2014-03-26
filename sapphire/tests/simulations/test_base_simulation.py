from mock import sentinel, Mock, patch, MagicMock
import unittest

import tables

from sapphire.simulations.base import BaseSimulation
from sapphire import storage


class BaseSimulationTest(unittest.TestCase):
    
    @patch.object(BaseSimulation, '_prepare_output_tables')
    def setUp(self, mock_method):
        self.mock_prepare_output_tables = mock_method
        self.cluster = sentinel.cluster
        self.datafile = sentinel.datafile
        self.output_path = sentinel.output_path
        self.N = sentinel.N

        self.simulation = BaseSimulation(self.cluster, self.datafile,
                                         self.output_path, self.N)

    def test_init_sets_attributes(self):
        self.assertIs(self.simulation.cluster, self.cluster)
        self.assertIs(self.simulation.datafile, self.datafile)
        self.assertIs(self.simulation.output_path, self.output_path)
        self.assertIs(self.simulation.N, self.N)

    def test_init_calls_prepare_output_tables(self):
        self.mock_prepare_output_tables.assert_called_once_with()

    @patch.object(BaseSimulation, '_prepare_coincidence_tables')
    @patch.object(BaseSimulation, '_prepare_station_tables')
    @patch.object(BaseSimulation, '_store_station_index')
    def test_prepare_output_tables_calls(self, mock_method3, mock_method2,
                                         mock_method1):
        self.simulation._prepare_output_tables()
        mock_method1.assert_called_once_with()
        mock_method2.assert_called_once_with()
        mock_method3.assert_called_once_with()

    @unittest.skip("WIP")
    def test_run(self):
        pass

    @unittest.skip("WIP")
    def test_generate_shower_parameters(self):
        pass

    @unittest.skip("WIP")
    def test_simulate_station_response(self, station, shower_parameters):
        pass

    @unittest.skip("WIP")
    def test_simulate_detector_response(self, detector, shower_parameters):
        pass

    @unittest.skip("WIP")
    def test_simulate_trigger(self, detector_observables):
        pass

    @unittest.skip("WIP")
    def test_simulate_gps(self, station_observables, shower_parameters, station):
        pass

    @unittest.skip("WIP")
    def test_process_detector_observables(self, detector_observables):
        pass

    @unittest.skip("WIP")
    def test_store_station_observables(self, station_id, station_observables):
        pass

    @unittest.skip("WIP")
    def test_store_coincidence(self, shower_id, shower_parameters, station_events):
        pass

    @unittest.skip("WIP")
    def test_prepare_coincidence_tables(self):
        pass

    @unittest.skip("WIP")
    def test_prepare_station_tables(self):
        pass

    @unittest.skip("WIP")
    def test_store_station_index(self):
        pass

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


@patch.object(BaseSimulation, 'generate_shower_parameters')
class BaseSimulationRunMethodTest(unittest.TestCase):

    @patch.object(BaseSimulation, '_prepare_output_tables')
    def setUp(self, mock_prepare_output_tables):
        self.simulation = BaseSimulation(Mock(), Mock(), Mock(), Mock())

    def test_run(self, mock_generate_shower_parameters):
        self.simulation.run()


if __name__ == '__main__':
    unittest.main()
