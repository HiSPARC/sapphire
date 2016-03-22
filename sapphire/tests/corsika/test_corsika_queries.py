import unittest
import os

from mock import sentinel, patch, MagicMock

from sapphire.corsika import corsika_queries

TEST_OVERVIEW_FILE = 'test_data/corsika_overview.h5'


class CorsikaQueryTest(unittest.TestCase):

    def setUp(self):
        self.cq = corsika_queries.CorsikaQuery(self.get_overview_path())

    def tearDown(self):
        self.cq.finish()

    def test_seeds(self):
        result = self.cq.seeds(self.cq.all_simulations())
        self.assertEqual(result, ['1_2'])

        result = self.cq.seeds(self.cq.all_simulations(), iterator=True)
        self.assertEqual(list(result), ['1_2'])

    def test_get_info(self):
        result = self.cq.get_info('1_2')
        self.assertEqual(result, self.cq.sims[0])
        self.assertRaises(ValueError, self.cq.get_info, '1')
        self.assertRaises(ValueError, self.cq.get_info, '1_2_3')
        self.assertRaises(IndexError, self.cq.get_info, '1_3')

    def test_all_energies(self):
        energies = self.cq.all_energies
        self.assertEqual(energies, set([14.]))

    def test_all_particles(self):
        particles = self.cq.all_particles
        self.assertEqual(particles, set(['proton']))

    def test_all_azimuths(self):
        azimuths = self.cq.all_azimuths
        self.assertEqual(azimuths, set([-90.]))

    def test_all_zeniths(self):
        zeniths = self.cq.all_zeniths
        self.assertEqual(zeniths, set([0.]))

    def test_available_parameters(self):
        result = self.cq.available_parameters('energy', particle='proton')
        self.assertEqual(result, {14.0})
        result = self.cq.available_parameters('particle_id', zenith=0.)
        self.assertEqual(result, {'proton'})
        result = self.cq.available_parameters('zenith', azimuth=-90.)
        self.assertEqual(result, {0.})
        self.assertRaises(RuntimeError, self.cq.available_parameters, 'zenith',
                          energy=19)
        self.assertRaises(RuntimeError, self.cq.available_parameters, 'zenith',
                          particle='iron')

    def get_overview_path(self):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, TEST_OVERVIEW_FILE)


class MockCorsikaQueryTest(unittest.TestCase):

    @patch.object(corsika_queries.tables, 'open_file')
    def setUp(self, mock_open):
        self.mock_open = mock_open
        self.data_path = sentinel.data_path
        self.simulations_group = sentinel.simulations_group

        self.cq = corsika_queries.CorsikaQuery(self.data_path,
                                               self.simulations_group)

    def test_init(self):
        self.mock_open.assert_called_once_with(self.data_path, 'r')
        self.mock_open.return_value.get_node.assert_called_once_with(
            sentinel.simulations_group)

    @patch.object(corsika_queries.tables, 'open_file')
    def test_init_file(self, mock_open):
        data = MagicMock(spec=corsika_queries.tables.File)
        corsika_queries.CorsikaQuery(data, sentinel.simulations_group)
        data.get_node.assert_called_once_with(sentinel.simulations_group)
        self.assertFalse(mock_open.called)

    def test_finish(self):
        self.cq.finish()
        self.cq.data.close.assert_called_once_with()

    @patch.object(corsika_queries.CorsikaQuery, 'perform_query')
    def test_simulations(self, mock_perform):
        mock_perform.return_value = sentinel.simulations
        result = self.cq.simulations(particle=None)
        self.assertEqual(result, sentinel.simulations)
        mock_perform.assert_called_once_with('', False)

        self.cq.all_particles = ['electron']
        self.cq.all_energies = [15.5]
        result = self.cq.simulations(particle='electron', energy=15.5,
                                     zenith=22.5, azimuth=90.)
        self.assertEqual(result, sentinel.simulations)
        mock_perform.assert_called_with(
            '(particle_id == 3) & '
            '(abs(log10(energy) - 15.5) < 1e-4) & '
            '(abs(zenith - 0.392699081699) < 1e-4) & '
            '(abs(azimuth - 1.57079632679) < 1e-4)', False)

    def test_filter(self):
        filter = self.cq.filter('type', 123)
        self.assertEqual(filter, '(type == 123)')

    def test_float_filter(self):
        filter = self.cq.float_filter('type', 12.3)
        self.assertEqual(filter, '(abs(type - 12.3) < 1e-4)')

    def test_range_filter(self):
        filter = self.cq.range_filter('type', 12.3, 14.5)
        self.assertEqual(filter, '(type >= 12.3) & (type <= 14.5)')
        filter = self.cq.range_filter('type', 12.3)
        self.assertEqual(filter, '(type >= 12.3)')
        filter = self.cq.range_filter('type', max=14.5)
        self.assertEqual(filter, '(type <= 14.5)')
        filter = self.cq.range_filter('type')
        self.assertEqual(filter, '')

    def test_all_simulations(self):
        result = self.cq.all_simulations()
        self.cq.sims.read.assert_called_once_with()
        self.assertEqual(result, self.cq.sims.read.return_value)

        result = self.cq.all_simulations(iterator=True)
        self.cq.sims.iterrows.assert_called_once_with()
        self.assertEqual(result, self.cq.sims.iterrows.return_value)

    def test_perform_query(self):
        result = self.cq.perform_query(sentinel.query, iterator=True)
        self.assertEqual(result, self.cq.sims.where.return_value)
        self.cq.sims.where.assert_called_once_with(sentinel.query)

        result = self.cq.perform_query(sentinel.query, iterator=False)
        self.assertEqual(result, self.cq.sims.read_where.return_value)
        self.cq.sims.read_where.assert_called_once_with(sentinel.query)


if __name__ == '__main__':
    unittest.main()
