import unittest
import os

from mock import sentinel, patch, call

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

    def test_finish(self):
        self.cq.finish()
        self.cq.data.close.assert_called_once_with()

    def test_filter(self):
        filter = self.cq.filter('type', 123)
        self.assertEqual(filter, '(type == 123)')

    def test_float_filter(self):
        filter = self.cq.float_filter('type', 12.3)
        self.assertEqual(filter, '(abs(type - 12.3) < 1e-5)')

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


if __name__ == '__main__':
    unittest.main()
