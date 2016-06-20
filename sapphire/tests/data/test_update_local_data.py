import unittest

from mock import patch

from sapphire.data import update_local_data


class UpdateLocalDataTests(unittest.TestCase):

    @patch.object(update_local_data, 'update_sublevel_json')
    @patch.object(update_local_data, 'update_toplevel_json')
    def test_update_local_json(self, mock_top, mock_sub):
        update_local_data.update_local_json(progress=False)
        self.assertTrue(mock_top.called)
        self.assertTrue(mock_sub.called)

    @patch.object(update_local_data, 'Network')
    @patch.object(update_local_data, 'HiSPARCNetwork')
    @patch.object(update_local_data, 'update_subsublevel_tsv')
    @patch.object(update_local_data, 'update_sublevel_tsv')
    def test_update_local_tsv(self, mock_sub, mock_ssub, mock_hnet, mock_net):
        update_local_data.update_local_tsv(progress=False)
        self.assertTrue(mock_sub.called)
        self.assertTrue(mock_ssub.called)


if __name__ == '__main__':
    unittest.main()
