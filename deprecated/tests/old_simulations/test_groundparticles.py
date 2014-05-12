import unittest
import types

import tables
from numpy import array, deg2rad
from math import pi, atan2, sqrt
from mock import Mock, MagicMock, patch, sentinel

from sapphire.old_simulations import groundparticles, BaseSimulation
from sapphire import clusters
from sapphire import storage


class GroundParticleSimulationTests(unittest.TestCase):
    @patch('os.path.split')
    def setUp(self, os_path_split_mock):
        self.output_head, self.output_tail = (sentinel.output_head,
                                              sentinel.output_tail)
        os_path_split_mock.return_value = self.output_head, self.output_tail

        self.cluster = sentinel.cluster
        self.data = MagicMock(name='data')
        self.grdpcles = Mock(name='grdpcles')
        self.grdpcles._v_pathname = '/E_1PeV/zenith_22_5/shower_0'
        self.data.get_node.return_value = self.grdpcles
        self.output = sentinel.output
        self.R = sentinel.R
        self.N = sentinel.N
        self.simulation = groundparticles.GroundParticlesSimulation(self.cluster, self.data,
                                              self.grdpcles._v_pathname, self.output,
                                              self.R, self.N)

        # make progressbar(list) do nothing (i.e., return list)
        self.progressbar_patcher = patch('progressbar.ProgressBar')
        self.progressbar_mock = self.progressbar_patcher.start()
        self.progressbar_mock.return_value.side_effect = lambda x: x

    def tearDown(self):
        self.progressbar_patcher.stop()

    def test_init_sets_attributes(self):
        self.assertIs(self.simulation.cluster, self.cluster)
        self.assertIs(self.simulation.data, self.data)
        self.assertIs(self.simulation.R, self.R)
        self.assertIs(self.simulation.N, self.N)

    def test_init_gets_grdpcles_node(self):
        self.data.get_node.assert_called_with(self.grdpcles._v_pathname)
        self.assertIs(self.simulation.grdpcles, self.data.get_node.return_value)

    @patch('os.path.split')
    def test_init_raises_runtimeerror_if_grdpcles_not_found(self, os_path_split_mock):
        os_path_split_mock.return_value = Mock(), Mock()

        data = MagicMock()
        data.get_node.side_effect = tables.NoSuchNodeError
        with self.assertRaises(RuntimeError):
            groundparticles.GroundParticlesSimulation(Mock(), data, Mock(), Mock(), Mock(), Mock())

    def test_init_creates_output_group(self):
        self.data.create_group.assert_called_with(self.output_head,
                                                 self.output_tail,
                                                 createparents=True)
        self.assertIs(self.simulation.output, self.data.create_group.return_value)

    def test_generate_positions_is_generator(self):
        self.assertEqual(type(self.simulation.generate_positions()), types.GeneratorType)

    def test_generate_positions_returns_uniform_circle_distribution(self):
        self.simulation.R = 100
        self.simulation.N = 10000
        positions = array(list(self.simulation.generate_positions()))
        r = positions[:, 0]
        phi = positions[:, 1]
        alpha = positions[:, 2]

        self.assertLessEqual(max(r), self.simulation.R)
        self.assertLessEqual(max(phi), pi)
        self.assertGreaterEqual(min(phi), -pi)
        self.assertLessEqual(max(alpha), pi)
        self.assertGreaterEqual(min(alpha), -pi)

        self.assertAlmostEqual(max(r), self.simulation.R, delta=1)
        self.assertAlmostEqual(max(phi), pi, delta=.01)
        self.assertAlmostEqual(min(phi), -pi, delta=.01)
        self.assertAlmostEqual(max(alpha), pi, delta=.01)
        self.assertAlmostEqual(min(alpha), -pi, delta=.01)

        # More distribution tests? X and Y, e.g.

    def test_get_station_particles(self):
        # mock station
        detectors = [sentinel.detector1, sentinel.detector2]
        station = Mock()
        station.detectors = detectors

        # mock get_detector_particles
        results = [sentinel.result2, sentinel.result1]
        get_detector_particles = Mock(side_effect=lambda * args: results.pop())
        self.simulation.get_detector_particles = get_detector_particles

        particles = self.simulation.get_station_particles(station)

        # assertions
        get_detector_particles.assert_called_with(detectors[1])
        pop_last_call(get_detector_particles)
        get_detector_particles.assert_called_with(detectors[0])
        self.assertEqual(particles, [sentinel.result1, sentinel.result2])

    def test_get_detector_particles(self):
        # I'm not terribly sure how much this helps...

        detector = Mock()
        corners = [Mock(), Mock(), Mock(), Mock()]
        detector.get_corners.return_value = corners

        get_line_boundary_eqs = Mock()
        eqs_results = [(0., 'y - x', 1.), (2., 'y + x', 3.)]
        get_line_boundary_eqs.side_effect = lambda * args: eqs_results.pop()
        self.simulation.get_line_boundary_eqs = get_line_boundary_eqs

        results = self.simulation.get_detector_particles(detector)
        self.simulation.grdpcles.read_where.assert_called_with(
            "(b11 < y + x) & (y + x < b12) & (b21 < y - x) & (y - x < b22)")
        self.assertIs(results, self.simulation.grdpcles.read_where.return_value)

        get_line_boundary_eqs.assert_called_with(corners[1], corners[2], corners[3])
        pop_last_call(get_line_boundary_eqs)
        get_line_boundary_eqs.assert_called_with(corners[0], corners[1], corners[2])

    def test_get_line_boundary_eqs(self):
        func = self.simulation.get_line_boundary_eqs
        self.assertEqual(func((0, 0), (1, 0), (1, 1)), (0, 'y - 0.000000 * x', 1))
        self.assertEqual(func((1, 0), (0, 0), (1, 1)), (0, 'y - 0.000000 * x', 1))
        self.assertEqual(func((1, 0), (0, 0), (1, -1)), (-1, 'y - 0.000000 * x', 0))

        self.assertEqual(func((0, 0), (1, 1), (0, 2)), (0, 'y - 1.000000 * x', 2))
        self.assertEqual(func((1, 1), (0, 0), (0, 2)), (0, 'y - 1.000000 * x', 2))

        self.assertEqual(func((0, 0), (0, 1), (2, 2)), (0, 'x', 2))
        self.assertEqual(func((0, 0), (0, 1), (-2, 2)), (-2, 'x', 0))

    def test_get_detector_particles_with_real_data(self):
        simulation = self.setup_simulation_with_real_data()

        station = Mock()
        detector = clusters.Detector(station, 0, 0, 'UD')

        station.get_xyalpha_coordinates.return_value = (0, 0, 0)
        particles = simulation.get_detector_particles(detector)
        ids = [x['id'] for x in particles]
        self.assertEqual(ids, [0, 1, 2, 5])

        station.get_xyalpha_coordinates.return_value = (.2, 0, 0)
        particles = simulation.get_detector_particles(detector)
        ids = [x['id'] for x in particles]
        self.assertEqual(ids, [0, 1, 3, 5])

        station.get_xyalpha_coordinates.return_value = (.2, 0, -pi / 4)
        particles = simulation.get_detector_particles(detector)
        ids = [x['id'] for x in particles]
        self.assertEqual(ids, [0, 3, 5])

        simulation.data.close()

    def setup_simulation_with_real_data(self):
        data = self.setup_datafile_with_real_data()
        return groundparticles.GroundParticlesSimulation(sentinel.cluster, data,
                                                         '/zenith_22_5/grdpcles',
                                                         '/output', sentinel.R, sentinel.N)

    def setup_datafile_with_real_data(self):
        data = tables.open_file('/tmp/tmp.h5', 'w')
        grdpcles = data.create_table('/zenith_22_5', 'grdpcles', storage.ShowerParticle,
                                    createparents=True)
        self.create_particle(grdpcles, 0, 0., 0.)
        self.create_particle(grdpcles, 1, .24, .49)
        self.create_particle(grdpcles, 2, -.24, -.49)
        self.create_particle(grdpcles, 3, .26, 0.)
        self.create_particle(grdpcles, 4, 0., .51)
        self.create_particle(grdpcles, 5, .2, .2)
        grdpcles.flush()
        return data

    def create_particle(self, table, id, x, y):
        row = table.row
        row['id'] = id
        row['x'] = x
        row['y'] = y
        row['core_distance'] = sqrt(x ** 2 + y ** 2)
        row['polar_angle'] = atan2(y, x)
        row.append()


    def test_run_without_arguments(self):
        my_tables = [Mock(), Mock()]
        self.data.create_table.side_effect = lambda * args: my_tables.pop()

        self.simulation._run_welcome_msg = Mock()
        self.simulation._run_exit_msg = Mock()
        self.simulation.generate_positions = Mock()
        self.simulation.generate_positions.return_value = [(0, 0, 0), (1, 1, 1)]
        self.simulation.simulate_event = Mock()
        self.simulation.store_observables = Mock()
        self.simulation.run()

        self.data.create_table.assert_called_with(self.simulation.output, '_particles', storage.SimulationParticle)
        pop_last_call(self.data.create_table)
        self.data.create_table.assert_called_with(self.simulation.output, '_headers', storage.SimulationEventHeader)
        self.simulation.headers.flush.assert_called_with()
        self.simulation.particles.flush.assert_called_with()
        self.simulation.store_observables.assert_called_with()

    def test_write_observables_uses_station_coordinates(self):
        observables = {'station_id': 1}
        t = [[1, 2], [2, 3, 4], [3, 4], [4, 5, 6]]
        num_particles = [len(u) for u in t]

        cluster = Mock(name='cluster')
        cluster.stations = MagicMock(name='stations')
        station = Mock(name='station')
        station.get_rphialpha_coordinates.return_value = sentinel.r, sentinel.phi, sentinel.alpha
        cluster.stations.__getitem__.return_value = station
        self.simulation.cluster = cluster

        with patch.object(BaseSimulation, 'write_observables') as mock_parent:
            self.simulation.write_observables(observables, num_particles, t)

        # list index must be zero-based
        cluster.stations.__getitem__.assert_called_once_with(0)
        # must get station coordinates
        station.get_rphialpha_coordinates.assert_called_once_with()
        # coordinates must be the station coordinates
        self.assertIs(observables['r'], sentinel.r)
        self.assertIs(observables['phi'], sentinel.phi)
        self.assertIs(observables['alpha'], sentinel.alpha)
        # must call parent method
        mock_parent.assert_called_once_with(observables, num_particles, t)

    def test_write_coincidence_calls_super_with_transformed_coordinates(self):
        event = {'phi': 1., 'alpha': 2.}
        transformed_event = {'phi': 1. + pi - 2., 'shower_theta': deg2rad(22.5),
                             'shower_phi': -2.}

        with patch.object(BaseSimulation, 'write_coincidence') as mock_parent:
            self.simulation.write_coincidence(event, sentinel.N)

        mock_parent.assert_called_once_with(transformed_event, sentinel.N)


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
