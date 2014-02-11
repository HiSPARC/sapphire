from mock import sentinel, Mock, patch, MagicMock
import unittest
from math import pi, sqrt

from sapphire.old_simulations.ldf import BaseLdfSimulation
from sapphire.old_simulations.base import BaseSimulation
import sapphire.clusters


class BaseLdfSimulationTest(unittest.TestCase):
    def setUp(self):
        cluster = Mock()
        data = MagicMock()
        output = '/simulations'
        R = sentinel.R
        N = 12345
        shower_size = sentinel.shower_size

        self.simulation = BaseLdfSimulation(cluster, data, output, R, N, shower_size)

        # make progressbar(list) do nothing (i.e., return list)
        self.progressbar_patcher = patch('progressbar.ProgressBar')
        self.progressbar_mock = self.progressbar_patcher.start()
        self.progressbar_mock.return_value.side_effect = lambda x: x

    def test_run_generates_positions(self):
        self.simulation._run_welcome_msg = Mock()
        self.simulation._run_exit_msg = Mock()
        self.simulation.observables.nrows = 0
        self.simulation.generate_positions = Mock()
        self.simulation.generate_positions.return_value = MagicMock()
        self.simulation.run()
        self.simulation.generate_positions.assert_called_with()

    def test_run_uses_positions(self):
        self.simulation._run_welcome_msg = Mock()
        self.simulation._run_exit_msg = Mock()
        self.simulation.observables.nrows = 0
        self.simulation.generate_positions = Mock()
        self.simulation.run(MagicMock())
        self.assertFalse(self.simulation.generate_positions.called)

    def test_run(self):
        self.simulation._run_welcome_msg = Mock()
        self.simulation._run_exit_msg = Mock()
        self.simulation.observables.nrows = 0
        self.simulation.simulate_event = Mock()

        self.simulation.run([(sentinel.r, sentinel.phi)])

        event = {'id': 0, 'r': sentinel.r, 'phi': sentinel.phi, 'alpha': 0.,
                 'shower_theta': 0., 'shower_phi': 0., 'shower_size': sentinel.shower_size}
        self.simulation.simulate_event.assert_called_with(event)

    def test_simulate_event(self):
        self.simulation.simulate_station_observables_and_return_has_triggered_and_eventid = Mock()
        self.simulation.simulate_station_observables_and_return_has_triggered_and_eventid.return_value = True, sentinel.station_event_id
        self.simulation.write_coincidence = Mock()

        self.simulation.cluster.stations = [sentinel.station]
        event = MagicMock()

        self.simulation.simulate_event(event)

        self.simulation.simulate_station_observables_and_return_has_triggered_and_eventid.assert_called_with(sentinel.station, event)
        self.simulation.write_coincidence.assert_called_with(event, 1)
        self.simulation.c_index.append.assert_called_with([sentinel.station_event_id])

    def test_simulate_event_calculates_correct_multiplicity(self):
        event = MagicMock()
        self.simulation.cluster.stations = MagicMock()
        self.simulation.write_coincidence = Mock()
        self.simulation.cluster.stations = [Mock(), Mock(), Mock()]
        has_triggered = SideEffects([(1, 1), (0, 2), (1, 3)])
        self.simulation.simulate_station_observables_and_return_has_triggered_and_eventid = Mock()
        self.simulation.simulate_station_observables_and_return_has_triggered_and_eventid.side_effect = has_triggered

        self.simulation.simulate_event(event)

        self.simulation.write_coincidence.assert_called_once_with(event, 2)

    def test_simulate_station_observables_and_return_has_triggered_and_eventid(self):
        self.simulation.simulate_detector_observables = Mock()
        num_particles = SideEffects([sentinel.n1, sentinel.n2, sentinel.n3, sentinel.n4])
        self.simulation.simulate_detector_observables.side_effect = num_particles
        self.simulation.write_observables_and_return_id = Mock()
        self.simulation.write_observables_and_return_id.return_value = 28

        station = Mock()
        station.detectors = [sentinel.detector1, sentinel.detector2,
                             sentinel.detector3, sentinel.detector4]
        event = sentinel.event

        has_triggered, observables_id = self.simulation.simulate_station_observables_and_return_has_triggered_and_eventid(station, event)

        self.assertEqual(has_triggered, True)
        self.assertEqual(observables_id, 28)

        self.simulation.simulate_detector_observables.assert_called_with(sentinel.detector4, sentinel.event)
        pop_last_call(self.simulation.simulate_detector_observables)
        self.simulation.simulate_detector_observables.assert_called_with(sentinel.detector3, sentinel.event)
        pop_last_call(self.simulation.simulate_detector_observables)
        self.simulation.simulate_detector_observables.assert_called_with(sentinel.detector2, sentinel.event)
        pop_last_call(self.simulation.simulate_detector_observables)
        self.simulation.simulate_detector_observables.assert_called_once_with(sentinel.detector1, sentinel.event)
        self.simulation.write_observables_and_return_id.assert_called_with(station, event, [sentinel.n1, sentinel.n2, sentinel.n3, sentinel.n4])

    @unittest.skip("Needs better mocking, or just drop it")
    def test_simulate_detector_observables(self):
        area = .5
        detector = Mock()
        detector.get_area.return_value = area

        N = 4.2
        self.simulation.calculate_core_distance = Mock()
        self.simulation.calculate_core_distance.return_value = sentinel.R
        self.simulation.calculate_ldf_value = Mock()
        self.simulation.calculate_ldf_value.return_value = N

        event = MagicMock()
        event.__getitem__.return_value = sentinel.shower_size

        num_particles = self.simulation.simulate_detector_observables(detector, event)

        self.assertEqual(num_particles, N * area)
        self.simulation.calculate_core_distance.assert_called_once_with(detector, event)
        self.simulation.calculate_ldf_value.assert_called_once_with(sentinel.R, sentinel.shower_size)
        detector.get_area.assert_called_once_with()
        event.__getitem__.assert_called_once_with('shower_size')

    @unittest.skip('Need better test')
    def test_calculate_core_distance(self):
        detector = Mock()
        detector.get_xy_coordinates.return_value = (4, 5)
        event = {'r': sqrt(2), 'phi': 3 * pi / 4}

        distance = sqrt((4 - -1) ** 2 + (5 - 1) ** 2)

        R = self.simulation.calculate_core_distance(detector, event)

        self.assertEqual(R, distance)
        detector.get_xy_coordinates.assert_called_once_with()

    def test_calculate_ldf_value(self):
        """The base class should NOT return particles"""

        value = self.simulation.calculate_ldf_value(sentinel.R, sentinel.shower_size)
        self.assertEqual(value, 0.)

    @unittest.skip("Broken test, FAIL")
    def test_write_observables_and_return_id(self):
        station = Mock()
        station.station_id = sentinel.station_id
        station.get_xyalpha_coordinates.return_value = sentinel.x, sentinel.y, sentinel.alpha
        station.get_rphialpha_coordinates.return_value = sentinel.r, sentinel.phi, sentinel.alpha
        event = MagicMock()
        event.__getitem__.return_value = sentinel.event_id
        n1, n2, n3, n4 = 5, 6, 0, 8

        row = self.simulation.observables.row
        self.simulation._observables_nrows = 27

        id = self.simulation.write_observables_and_return_id(station, event, [n1, n2, n3, n4])

        self.assertEqual(id, 27)
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'id', sentinel.event_id))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'station_id', sentinel.station_id))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'r', sentinel.r))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'phi', sentinel.phi))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'x', sentinel.x))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'y', sentinel.y))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'alpha', sentinel.alpha))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'N', 3))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 't1', 0))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 't2', 0))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 't3', 0))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 't4', 0))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'n1', 5))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'n2', 6))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'n3', 0))
        self.assertTrue(is_mock_previously_called_with(row.__setitem__, 'n4', 8))
        self.simulation.observables.row.append.assert_called_once_with()
        self.assertEqual(self.simulation._observables_nrows, 28)


class SideEffects:
    def __init__(self, response_list):
        self.response_list = response_list

    def __call__(self, *args):
        return self.response_list.pop(0)


def is_mock_previously_called_with(mock, *args, **kwargs):
    return (args, kwargs) in mock.call_args_list


def pop_last_call(mock):
    if not mock.call_count:
        raise AssertionError("Cannot pop last call: call_count is 0")
    mock.call_args_list.pop()
    try:
        mock.call_args = mock.call_args_list[-1]
    except IndexError:
        mock.call_args = None
        mock.called = False
    mock.call_count -= 1


if __name__ == '__main__':
    unittest.main()
