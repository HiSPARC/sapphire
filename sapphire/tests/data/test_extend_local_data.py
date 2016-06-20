import unittest

from mock import patch, sentinel

from sapphire.data import extend_local_data


class UpdateLocalDataTests(unittest.TestCase):

    @patch.object(extend_local_data, 'Network')
    @patch.object(extend_local_data, 'update_sublevel_tsv')
    def test_update_local_json(self, mock_sub, mock_net):
        mock_net.return_value.station_numbers.return_value = sentinel.numbers
        extend_local_data.update_additional_local_tsv(progress=False)
        mock_sub.assert_called_once_with('eventtime', sentinel.numbers, False)


if __name__ == '__main__':
    unittest.main()
