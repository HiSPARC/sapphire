import unittest
import warnings
import random

from mock import patch, sentinel, mock_open

from sapphire.corsika import qsub_corsika


class CorsikaBatchTest(unittest.TestCase):

    def setUp(self):
        self.cb = qsub_corsika.CorsikaBatch()

    @patch.object(qsub_corsika.particles, 'particle_id')
    def test_init(self, mock_particles):
        mock_particles.return_value = sentinel.particle_id
        cb = qsub_corsika.CorsikaBatch(16, sentinel.particle, sentinel.zenith,
                                       30, sentinel.queue, sentinel.corsika)
        self.assertEqual(cb.energy_pre, 1.)
        self.assertEqual(cb.energy_pow, 7)
        mock_particles.assert_called_once_with(sentinel.particle)
        self.assertEqual(cb.particle, sentinel.particle_id)
        self.assertEqual(cb.theta, sentinel.zenith)
        self.assertEqual(cb.phi, 120)
        self.assertEqual(cb.queue, sentinel.queue)
        self.assertEqual(cb.corsika, sentinel.corsika)
        self.assertIsNone(cb.seed1)
        self.assertIsNone(cb.seed2)
        self.assertIsNone(cb.rundir)

    def test_init_fractional_energy(self):
        cb = qsub_corsika.CorsikaBatch(16.5)
        self.assertEqual(cb.energy_pre, 3.16228)
        self.assertEqual(cb.energy_pow, 7)

    def test_init_bad_energy(self):
        self.assertRaises(ValueError, qsub_corsika.CorsikaBatch, 16.15)

    @patch.object(qsub_corsika.CorsikaBatch, 'prepare_env')
    @patch.object(qsub_corsika.CorsikaBatch, 'submit_job')
    def test_run(self, mock_submit_job, mock_prepare_env):
        self.cb.run()
        mock_prepare_env.assert_called_once_with()
        mock_submit_job.assert_called_once_with()

    @patch.object(qsub_corsika.qsub, 'submit_job')
    @patch.object(qsub_corsika.CorsikaBatch, 'get_rundir')
    @patch.object(qsub_corsika.CorsikaBatch, 'create_script')
    def test_submit_job(self, mock_create_script, mock_rundir, mock_submit_job):
        self.cb.seed1 = 123
        self.cb.seed2 = 456
        mock_rundir.return_value = '/data/123_456/'
        mock_create_script.return_value = sentinel.script
        self.cb.submit_job()
        mock_submit_job.assert_called_once_with(sentinel.script, 'cor_123_456',
                                                'generic', '-d /data/123_456/')
        # Check addition of walltime argument for long queue
        self.cb.queue = 'long'
        self.cb.submit_job()
        mock_submit_job.assert_called_with(sentinel.script, 'cor_123_456',
                                           'long',
                                           '-d /data/123_456/ -l walltime=96:00:00')

    @patch.object(qsub_corsika.os, 'listdir')
    def test_taken_seeds(self, mock_listdir):
        mock_listdir.return_value = [sentinel.dirs]
        taken = self.cb.taken_seeds()
        self.assertEqual(taken, [sentinel.dirs, sentinel.dirs])
        mock_listdir.assert_any_call(qsub_corsika.DATADIR)
        mock_listdir.assert_called_with(qsub_corsika.TEMPDIR)

    def test_generate_random_seeds(self):
        random.seed(0)
        self.cb.generate_random_seeds(['68764531_6546560', '32716_687164'])
        self.assertEqual(self.cb.seed1, 759979667)
        self.assertEqual(self.cb.seed2, 682158963)
        self.assertEqual(self.cb.rundir, '759979667_682158963/')

    def test_generated_random_seeds_taken(self):
        random.seed(0)
        self.cb.generate_random_seeds(['759979667_682158963', '32716_687164'])
        self.assertEqual(self.cb.seed1, 378514423)
        self.assertEqual(self.cb.seed2, 233025076)
        self.assertEqual(self.cb.rundir, '378514423_233025076/')

    @patch.object(qsub_corsika.os, 'mkdir')
    @patch.object(qsub_corsika.CorsikaBatch, 'get_rundir')
    def test_make_rundir(self, mock_rundir, mock_mkdir):
        mock_rundir.return_value = sentinel.rundir
        self.cb.make_rundir()
        mock_mkdir.assert_called_once_with(sentinel.rundir)

    @patch.object(qsub_corsika.os, 'chdir')
    @patch.object(qsub_corsika.CorsikaBatch, 'get_rundir')
    def test_goto_rundir(self, mock_rundir, mock_chdir):
        mock_rundir.return_value = sentinel.rundir
        self.cb.goto_rundir()
        mock_chdir.assert_called_once_with(sentinel.rundir)

    def test_get_rundir(self):
        self.cb.rundir = '123_456/'
        rundir = self.cb.get_rundir()
        self.assertEqual(rundir, qsub_corsika.TEMPDIR + '123_456/')

    @patch.object(qsub_corsika.CorsikaBatch, 'get_rundir')
    def test_create_input(self, mock_rundir):
        mock_rundir.return_value = '/data/123_456'
        mock_file = mock_open()
        with patch('__builtin__.open', mock_file):
            self.cb.create_input()
        mock_rundir.assert_called_once_with()
        mock_file.assert_called_once_with('/data/123_456/input-hisparc', 'w')
        self.assertTrue(mock_file().write.called)


class MultipleJobsTest(unittest.TestCase):

    @patch.object(qsub_corsika.qsub, 'check_queue')
    def test_no_available_slots(self, mock_check_queue):
        """No slots available on queue"""

        mock_check_queue.return_value = 0
        self.assertRaises(Exception, qsub_corsika.multiple_jobs, sentinel.n,
                          sentinel.energy, sentinel.particle, sentinel.zenith,
                          sentinel.azimuth, sentinel.queue, sentinel.corsika,
                          progress=False)
        mock_check_queue.assert_called_once_with(sentinel.queue)

    @patch.object(qsub_corsika, 'CorsikaBatch')
    @patch.object(qsub_corsika.qsub, 'check_queue')
    def test_one_available_wanted_more(self, mock_check_queue, mock_corsika_batch):
        """Only one slot available on queue"""

        mock_check_queue.return_value = 1
        with warnings.catch_warnings(record=True) as warned:
            qsub_corsika.multiple_jobs(2, sentinel.energy, sentinel.particle,
                                       sentinel.zenith, sentinel.azimuth,
                                       sentinel.queue, sentinel.corsika,
                                       progress=False)
        mock_check_queue.assert_called_once_with(sentinel.queue)
        mock_corsika_batch.assert_called_once_with(
            energy=sentinel.energy, particle=sentinel.particle,
            zenith=sentinel.zenith, azimuth=sentinel.azimuth,
            queue=sentinel.queue, corsika=sentinel.corsika)
        mock_corsika_batch.return_value.run.assert_called_once_with()
        self.assertEqual(len(warned), 1)

    @patch.object(qsub_corsika, 'CorsikaBatch')
    @patch.object(qsub_corsika.qsub, 'check_queue')
    def test_two_available_wanted_more(self, mock_check_queue, mock_corsika_batch):
        """Only two slots available on queue"""

        mock_check_queue.return_value = 2
        with warnings.catch_warnings(record=True) as warned:
            qsub_corsika.multiple_jobs(3, sentinel.energy, sentinel.particle,
                                       sentinel.zenith, sentinel.azimuth,
                                       sentinel.queue, sentinel.corsika,
                                       progress=False)
        mock_check_queue.assert_called_once_with(sentinel.queue)
        mock_corsika_batch.assert_called_with(
            energy=sentinel.energy, particle=sentinel.particle,
            zenith=sentinel.zenith, azimuth=sentinel.azimuth,
            queue=sentinel.queue, corsika=sentinel.corsika)
        mock_corsika_batch.return_value.run.assert_called_with()
        # This is twice as often because it includes the calls to run()
        self.assertEqual(len(mock_corsika_batch.mock_calls), 4)
        self.assertEqual(len(mock_corsika_batch.return_value.run.mock_calls), 2)
        self.assertEqual(len(warned), 1)

    @patch.object(qsub_corsika, 'CorsikaBatch')
    @patch.object(qsub_corsika.qsub, 'check_queue')
    def test_plenty_available(self, mock_check_queue, mock_corsika_batch):
        """Plenty of space on queue"""

        mock_check_queue.return_value = 50
        n = 10
        with warnings.catch_warnings(record=True) as warned:
            qsub_corsika.multiple_jobs(n, sentinel.energy, sentinel.particle,
                                       sentinel.zenith, sentinel.azimuth,
                                       sentinel.queue, sentinel.corsika,
                                       progress=False)
        mock_check_queue.assert_called_once_with(sentinel.queue)
        mock_corsika_batch.assert_called_with(
            energy=sentinel.energy, particle=sentinel.particle,
            zenith=sentinel.zenith, azimuth=sentinel.azimuth,
            queue=sentinel.queue, corsika=sentinel.corsika)
        mock_corsika_batch.return_value.run.assert_called_with()
        # This is twice as often because it includes the calls to run()
        self.assertEqual(len(mock_corsika_batch.mock_calls), n * 2)
        self.assertEqual(len(mock_corsika_batch.return_value.run.mock_calls), n)
        self.assertEqual(len(warned), 0)


if __name__ == '__main__':
    unittest.main()
