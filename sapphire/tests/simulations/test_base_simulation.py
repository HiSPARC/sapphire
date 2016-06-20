from mock import sentinel, Mock, patch, MagicMock, call
import types
import unittest
import warnings

import tables

from sapphire.simulations.base import BaseSimulation
from sapphire import storage


class BaseSimulationTest(unittest.TestCase):

    @patch.object(BaseSimulation, '_prepare_output_tables')
    def setUp(self, mock_method):
        self.mock_prepare_output_tables = mock_method
        self.cluster = sentinel.cluster
        self.data = sentinel.data
        self.output_path = sentinel.output_path
        self.N = sentinel.N

        self.simulation = BaseSimulation(self.cluster, self.data,
                                         self.output_path, self.N,
                                         progress=False)

    def test_init_sets_attributes(self):
        self.assertIs(self.simulation.cluster, self.cluster)
        self.assertIs(self.simulation.data, self.data)
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

    @patch.object(BaseSimulation, 'generate_shower_parameters')
    @patch.object(BaseSimulation, 'simulate_events_for_shower')
    @patch.object(BaseSimulation, 'store_coincidence')
    def test_run(self, mock_store, mock_simulate, mock_generate):
        mock_generate.return_value = [sentinel.params1, sentinel.params2]
        mock_simulate.return_value = sentinel.events
        self.simulation.run()

        # test simulate_events_for_shower called two times with
        # shower_parameters
        expected = [call(sentinel.params1), call(sentinel.params2)]
        self.assertEqual(mock_simulate.call_args_list, expected)

        # test store_coincidence called 2nd time with shower_id 1,
        # parameters and events
        mock_store.assert_called_with(1, sentinel.params2,
                                      sentinel.events)

    def test_generate_shower_parameters(self):
        self.simulation.N = 10
        output = self.simulation.generate_shower_parameters()
        self.assertIsInstance(output, types.GeneratorType)

        output = list(output)
        self.assertEqual(len(output), 10)

        expected = {'core_pos': (None, None), 'zenith': None, 'azimuth': None,
                    'size': None, 'energy': None, 'ext_timestamp': None}
        self.assertEqual(output[0], expected)

    @patch.object(BaseSimulation, 'simulate_station_response')
    @patch.object(BaseSimulation, 'store_station_observables')
    def test_simulate_events_for_shower(self, mock_store, mock_simulate):
        self.simulation.cluster = Mock()
        self.simulation.cluster.stations = [sentinel.station1,
                                            sentinel.station2,
                                            sentinel.station3]

        mock_simulate.side_effect = [(True, sentinel.obs1), (False, None),
                                     (True, sentinel.obs3)]
        mock_store.side_effect = [sentinel.index1, sentinel.index2]
        events = self.simulation.simulate_events_for_shower(
            sentinel.params)

        # test simulate_station_response called for each station, with
        # shower parameters
        expected = [call(sentinel.station1, sentinel.params),
                    call(sentinel.station2, sentinel.params),
                    call(sentinel.station3, sentinel.params)]
        self.assertEqual(mock_simulate.call_args_list, expected)

        # test store_station_observables called only for triggered
        # stations, with observables
        expected = [call(0, sentinel.obs1), call(2, sentinel.obs3)]
        self.assertEqual(mock_store.call_args_list, expected)

        # test returned events consists of list of station indexes and
        # stored event indexes
        self.assertEqual(events, [(0, sentinel.index1),
                                  (2, sentinel.index2)])

    @patch.object(BaseSimulation, 'simulate_all_detectors')
    @patch.object(BaseSimulation, 'simulate_trigger')
    @patch.object(BaseSimulation, 'process_detector_observables')
    @patch.object(BaseSimulation, 'simulate_gps')
    def test_simulate_station_response(self, mock_gps, mock_process,
                                       mock_trigger, mock_detectors):
        mock_detectors.return_value = sentinel.detector_observables
        mock_trigger.return_value = sentinel.has_triggered
        mock_process.return_value = sentinel.station_observables
        mock_gps.return_value = sentinel.gps_observables

        mock_station = Mock()
        mock_station.detectors = sentinel.detectors

        has_triggered, station_observables = \
            self.simulation.simulate_station_response(mock_station,
                                                      sentinel.parameters)

        # Tests
        mock_detectors.assert_called_once_with(sentinel.detectors,
                                               sentinel.parameters)
        mock_trigger.assert_called_once_with(sentinel.detector_observables)
        mock_process.assert_called_once_with(sentinel.detector_observables)
        mock_gps.assert_called_once_with(sentinel.station_observables,
                                         sentinel.parameters,
                                         mock_station)
        self.assertIs(has_triggered, sentinel.has_triggered)
        self.assertIs(station_observables, sentinel.gps_observables)

    @patch.object(BaseSimulation, 'simulate_detector_response')
    def test_simulate_all_detectors(self, mock_response):
        detectors = [sentinel.detector1, sentinel.detector2]
        mock_response.side_effect = [sentinel.observables1,
                                     sentinel.observables2]

        observables = self.simulation.simulate_all_detectors(
            detectors, sentinel.parameters)

        expected = [call(sentinel.detector1, sentinel.parameters),
                    call(sentinel.detector2, sentinel.parameters)]
        self.assertEqual(mock_response.call_args_list, expected)

        self.assertEqual(observables, [sentinel.observables1,
                                       sentinel.observables2])

    def test_simulate_detector_response(self):
        observables = self.simulation.simulate_detector_response(Mock(),
                                                                 Mock())
        self.assertIsInstance(observables, dict)
        self.assertIn('n', observables)
        self.assertIn('t', observables)

    def test_simulate_trigger(self):
        has_triggered = self.simulation.simulate_trigger(Mock())
        self.assertIsInstance(has_triggered, bool)

    def test_simulate_gps(self):
        mock_observables = Mock()
        self.simulation.simulate_gps(mock_observables, Mock(), Mock())
        self.assertEqual(mock_observables.update.call_count, 1)
        args, kwargs = mock_observables.update.call_args
        gps_dict = args[0]
        self.assertIsInstance(gps_dict, dict)
        self.assertIn('ext_timestamp', gps_dict)
        self.assertIn('timestamp', gps_dict)
        self.assertIn('nanoseconds', gps_dict)

    def test_process_detector_observables(self):
        detector_observables = [{'n': 1., 't': 2., 'pulseheights': 3.,
                                 'integrals': 4.},
                                {'n': 5., 't': 6., 'pulseheights': 7.,
                                 'integrals': 8.},
                                {'foo': -999.}]

        expected = {'n1': 1., 'n2': 5., 't1': 2., 't2': 6.,
                    'pulseheights': [3., 7., -1., -1.],
                    'integrals': [4., 8., -1, -1]}
        actual = self.simulation.process_detector_observables(
            detector_observables)

        self.assertEqual(expected, actual)

    def test_store_station_observables(self):
        station_groups = MagicMock()
        self.simulation.station_groups = station_groups
        table = station_groups.__getitem__.return_value.events
        table.nrows = 123

        observables = {'key1': 1., 'key2': 2.}
        table.colnames = ['key1', 'key2']
        idx = self.simulation.store_station_observables(
            sentinel.station_id, observables)

        # tests
        station_groups.__getitem__.assert_called_once_with(sentinel.station_id)

        expected = [call('event_id', table.nrows), call('key2', 2.),
                    call('key1', 1.)]
        self.assertEqual(table.row.__setitem__.call_args_list, expected)
        table.row.append.assert_called_once_with()
        table.flush.assert_called_once_with()
        self.assertEqual(idx, table.nrows - 1)

    def test_store_station_observables_raises_warning(self):
        station_groups = MagicMock()
        self.simulation.station_groups = station_groups
        table = station_groups.__getitem__.return_value.events
        observables = {'key1': 1., 'key2': 2.}
        table.colnames = ['key1']

        warnings.simplefilter('error')
        self.assertRaises(UserWarning,
                          self.simulation.store_station_observables,
                          sentinel.station_id, observables)

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
        self.data.create_group.assert_any_call(
            self.output_path, 'coincidences', createparents=True)
        self.data.create_table.assert_called_with(
            self.simulation.coincidence_group, 'coincidences', storage.Coincidence)
        self.assertEqual(self.data.create_vlarray.call_count, 2)
        self.data.create_vlarray.assert_any_call(
            self.simulation.coincidence_group, 'c_index', tables.UInt32Col(shape=2))

    @unittest.skip("Does not test this unit")
    def test_init_creates_cluster_output_group(self):
        self.data.create_group.assert_any_call(
            self.output_path, 'cluster_simulations', createparents=True)
        # The following tests need a better mock of cluster in order to work.
        # self.data.create_group.assert_any_call(self.simulation.cluster_group, 'station_0')
        # self.data.create_table.assert_any_call(
        #     station_group, 'events', storage.ProcessedHisparcEvent, expectedrows=self.N)

    @unittest.skip("Does not test this unit")
    def test_init_stores_cluster_in_attrs(self):
        self.assertIs(self.simulation.coincidence_group._v_attrs.cluster, self.cluster)


if __name__ == '__main__':
    unittest.main()
