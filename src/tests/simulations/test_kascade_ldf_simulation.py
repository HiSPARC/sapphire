import unittest
from mock import Mock, MagicMock, sentinel, patch
from numpy import logspace, array

from simulations.ldf import KascadeLdfSimulation, BaseLdfSimulation

class BaseLdfSimulationTest(unittest.TestCase):
    def setUp(self):
        cluster = Mock()
        data = MagicMock()
        output = '/simulations'
        R = sentinel.R
        N = MagicMock()

        self.simulation = KascadeLdfSimulation(cluster, data, output, R, N)

        # make progressbar(list) do nothing (i.e., return list)
        self.progressbar_patcher = patch('progressbar.ProgressBar')
        self.progressbar_mock = self.progressbar_patcher.start()
        self.progressbar_mock.return_value.side_effect = lambda x: x

    def test_cache_c_s_value(self):
        # change s value
        self.simulation._s += .1

        self.simulation.cache_c_s_value()
        expected = self.simulation._c(self.simulation._s)

        self.assertEqual(expected, self.simulation._c_s)

    def test_init_calls_cache_c_s_value(self):
        with patch.object(KascadeLdfSimulation, 'cache_c_s_value') as mock_cache:
            sim = KascadeLdfSimulation(Mock(), MagicMock(), '/mock', Mock(), Mock())
            mock_cache.assert_called_once_with()

    def test_init_calls_super_init_with_all_args(self):
        with patch.object(BaseLdfSimulation, '__init__') as mock_base_init:
            cluster = sentinel.cluster
            data = sentinel.data
            output = sentinel.output
            R = sentinel.R
            N = sentinel.N
            sim = KascadeLdfSimulation(cluster, data, output, R, N=N)
            mock_base_init.assert_called_once_with(cluster, data, output, R, N=N)

    def test_calculate_ldf_value_useful_immediately_after_init(self):
        self.simulation.calculate_ldf_value(1.)

    def test_run_updates_cached_c_s_value(self):
        sim = self.simulation

        sim._s += .1
        expected = sim._c(self.simulation._s)

        sim._run_welcome_msg = Mock()
        sim._run_exit_msg = Mock()
        sim.generate_positions = Mock()
        sim.generate_positions.return_value = MagicMock()
        sim.simulate_event = Mock()
        sim.observables.nrows = sim._observables_nrows

        sim.run()

        self.assertEqual(expected, sim._c_s)

    def test_run_calls_super_run(self):
        with patch.object(BaseLdfSimulation, 'run') as mock_run:
            self.simulation.run()
            mock_run.assert_called_once_with()

    def test_calculate_ldf_value(self):
        r = logspace(1, 3, 10)
        expected = array([  1.39434005e+01, 7.49964318e+00, 3.49327708e+00,
                            1.37707188e+00, 4.57377254e-01, 1.30067636e-01,
                            3.26255672e-02, 7.46150200e-03, 1.60164385e-03,
                            3.29852843e-04])
        actual = self.simulation.calculate_ldf_value(r)

        self.assertTrue(((expected - actual) / expected < 1e-8).all())
