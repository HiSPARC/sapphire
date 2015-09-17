import unittest
from itertools import cycle

from mock import patch, sentinel, mock_open

from sapphire import qsub


@patch.object(qsub.utils, 'which')
class CheckQueueTest(unittest.TestCase):

    @patch.object(qsub.subprocess, 'check_output')
    def test_queues(self, mock_check_output, mock_which):
        for queue in ['express', 'short', 'generic', 'long']:
            qsub.check_queue(queue)
            last_two_calls = mock_check_output.call_args_list[-2:]
            for call in last_two_calls:
                self.assertTrue(queue in call[0][0])

    @patch.object(qsub.subprocess, 'check_output')
    def test_bad_queue(self, mock_check_output, mock_which):
        self.assertRaises(KeyError, qsub.check_queue, 'bla')

    @patch.object(qsub.subprocess, 'check_output')
    def test_check_queue(self, mock_check_output, mock_which):
        combinations = ([['   0\n'], 2, 'express'],
                        [['   2\n'], 0, 'express'],
                        [[' 100\n'], 900, 'short'],
                        [['1100\n'], -100, 'short'],
                        [['2000\n', '1000\n'], 1000, 'generic'],
                        [['3600\n', '1000\n'], 400, 'generic'],
                        [[' 200\n', ' 100\n'], 400, 'long'],
                        [[' 620\n', ' 100\n'], 380, 'long'])
        for taken, available, queue in combinations:
            mock_check_output.side_effect = cycle(taken)
            self.assertEqual(qsub.check_queue(queue), available)


@patch.object(qsub.utils, 'which')
class SubmitJobTest(unittest.TestCase):

    @patch.object(qsub, 'create_script')
    @patch.object(qsub.subprocess, 'check_output')
    @patch.object(qsub, 'delete_script')
    def test_submit_job(self, mock_delete, mock_check_output, mock_create,
                        mock_which):
        mock_create.return_value = (sentinel.script_path, sentinel.script_name)
        mock_check_output.return_value = ''
        qsub.submit_job(sentinel.script, sentinel.name, sentinel.queue, sentinel.extra)

        mock_create.assert_called_once_with(sentinel.script, sentinel.name)
        command = ('qsub -q {queue} -V -z -j oe -N {name} {extra} {script}'
                   .format(queue=sentinel.queue, name=sentinel.script_name,
                           script=sentinel.script_path, extra=sentinel.extra))
        mock_check_output.assert_called_once_with(command,
                                                  stderr=qsub.subprocess.STDOUT,
                                                  shell=True)
        mock_delete.assert_called_once_with(sentinel.script_path)

    @patch.object(qsub, 'create_script')
    @patch.object(qsub.subprocess, 'check_output')
    @patch.object(qsub, 'delete_script')
    def test_failed_submit_job(self, mock_delete, mock_check_output,
                               mock_create, mock_which):
        mock_create.return_value = (sentinel.script_path, sentinel.script_name)
        mock_check_output.return_value = 'Failed!'
        self.assertRaises(Exception, qsub.submit_job, sentinel.script,
                          sentinel.name, sentinel.queue, sentinel.extra)

        mock_create.assert_called_once_with(sentinel.script, sentinel.name)
        command = ('qsub -q {queue} -V -z -j oe -N {name} {extra} {script}'
                   .format(queue=sentinel.queue, name=sentinel.script_name,
                           script=sentinel.script_path, extra=sentinel.extra))
        mock_check_output.assert_called_once_with(command,
                                                  stderr=qsub.subprocess.STDOUT,
                                                  shell=True)
        self.assertFalse(mock_delete.called)


class CreateScriptTest(unittest.TestCase):

    @patch.object(qsub.os, 'chmod')
    def test_create_script(self, mock_chmod):
        mock_file = mock_open()
        with patch('__builtin__.open', mock_file):
            res_path, res_name = qsub.create_script(sentinel.script, 'hoi')
        self.assertEqual(res_path, '/tmp/his_hoi.sh')
        self.assertEqual(res_name, 'his_hoi.sh')
        mock_file.assert_called_once_with(res_path, 'w')
        mock_file().write.called_once_with(sentinel.script)
        mock_chmod.assert_called_once_with(res_path, 0774)


class DeleteScriptTest(unittest.TestCase):

    @patch.object(qsub.os, 'remove')
    def test_delete_script(self, mock_remove):
        self.assertFalse(mock_remove.called)
        qsub.delete_script(sentinel.path)
        mock_remove.assert_called_once_with(sentinel.path)


if __name__ == '__main__':
    unittest.main()
