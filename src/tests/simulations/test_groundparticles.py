import unittest
import types

import tables
from numpy import array
from math import pi, atan2, sqrt
from mock import Mock, patch, sentinel

from simulations import groundparticles
import clusters
import storage

class BaseSimulationTests(unittest.TestCase):
    @patch('os.path.split')
    def setUp(self, os_path_split_mock):
        self.output_head, self.output_tail = (Mock(name='output_head'),
                                              Mock(name='output_tail'))
        os_path_split_mock.return_value = self.output_head, self.output_tail

        self.cluster = sentinel.cluster
        self.data = Mock(name='data')
        self.grdpcles = sentinel.grdpcles
        self.output = sentinel.output
        self.R = sentinel.R
        self.N = sentinel.N
        self.simulation = groundparticles.GroundParticlesSimulation(self.cluster, self.data,
                                              self.grdpcles, self.output,
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
        self.data.getNode.assert_called_with('/', self.grdpcles)
        self.assertIs(self.simulation.grdpcles, self.data.getNode.return_value)

    @patch('os.path.split')
    def test_init_raises_runtimeerror_if_grdpcles_not_found(self, os_path_split_mock):
        os_path_split_mock.return_value = Mock(), Mock()

        data = Mock()
        data.getNode.side_effect = tables.NoSuchNodeError
        with self.assertRaises(RuntimeError):
            groundparticles.GroundParticlesSimulation(Mock(), data, Mock(), Mock(), Mock(), Mock())

    def test_init_creates_output_group(self):
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
            groundparticles.GroundParticlesSimulation(Mock(), data, Mock(), Mock(), Mock(), Mock())

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

        # mock simple args
        X, Y = sentinel.X, sentinel.Y
        alpha = sentinel.alpha

        # mock get_detector_particles
        results = [sentinel.result2, sentinel.result1]
        get_detector_particles = Mock(side_effect=lambda * args: results.pop())
        self.simulation.get_detector_particles = get_detector_particles

        particles = self.simulation.get_station_particles(station, X, Y, alpha)

        # assertions
        get_detector_particles.assert_called_with(X, Y, detectors[1], alpha)
        pop_last_call(get_detector_particles)
        get_detector_particles.assert_called_with(X, Y, detectors[0], alpha)
        self.assertEqual(particles, [sentinel.result1, sentinel.result2])

    def test_get_detector_particles(self):
        # I'm not terribly sure how much this helps...

        X, Y, alpha = sentinel.X, sentinel.Y, sentinel.alpha
        detector = Mock()
        corners = [Mock(), Mock(), Mock(), Mock()]
        detector.get_corners.return_value = corners

        get_line_boundary_eqs = Mock()
        eqs_results = [(0., 'y - x', 1.), (2., 'y + x', 3.)]
        get_line_boundary_eqs.side_effect = lambda * args: eqs_results.pop()
        self.simulation.get_line_boundary_eqs = get_line_boundary_eqs

        results = self.simulation.get_detector_particles(X, Y, detector, alpha)
        self.simulation.grdpcles.readWhere.assert_called_with(
            "(b11 < y + x) & (y + x < b12) & (b21 < y - x) & (y - x < b22)")
        self.assertIs(results, self.simulation.grdpcles.readWhere.return_value)

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
        detector = clusters.Detector(sentinel.station, 0, 0, 'UD')

        particles = simulation.get_detector_particles(0, 0, detector, 0)
        ids = [x['id'] for x in particles]
        self.assertEqual(ids, [0, 1, 2, 5])

        particles = simulation.get_detector_particles(.2, 0, detector, 0)
        ids = [x['id'] for x in particles]
        self.assertEqual(ids, [0, 1, 3, 5])

        particles = simulation.get_detector_particles(.2, 0, detector, -pi / 4)
        ids = [x['id'] for x in particles]
        self.assertEqual(ids, [0, 3, 5])

        simulation.data.close()

    def setup_simulation_with_real_data(self):
        data = self.setup_datafile_with_real_data()
        return groundparticles.GroundParticlesSimulation(sentinel.cluster, data, '/grdpcles',
                                   '/output', sentinel.R, sentinel.N)

    def setup_datafile_with_real_data(self):
        data = tables.openFile('/tmp/tmp.h5', 'w')
        grdpcles = data.createTable('/', 'grdpcles', storage.Particle)
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
        self.data.createTable.side_effect = lambda * args: my_tables.pop()

        self.simulation._run_welcome_msg = Mock()
        self.simulation.generate_positions = Mock()
        self.simulation.generate_positions.return_value = [(0, 0, 0), (1, 1, 1)]
        self.simulation.simulate_event = Mock()
        self.simulation.store_observables = Mock()
        self.simulation.run()

        self.data.createTable.assert_called_with(self.simulation.output, '_particles', storage.ParticleEvent)
        pop_last_call(self.data.createTable)
        self.data.createTable.assert_called_with(self.simulation.output, '_headers', storage.SimulationHeader)
        self.simulation.headers.flush.assert_called_with()
        self.simulation.particles.flush.assert_called_with()
        self.simulation.store_observables.assert_called_with()


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
