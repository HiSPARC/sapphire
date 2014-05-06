from mock import sentinel, Mock, patch, MagicMock, call
import types
import unittest
import warnings

from sapphire.analysis import coincidence_queries


class BaseCoincidenceQueryTest(unittest.TestCase):

    @patch.object(coincidence_queries.tables, 'open_file')
    def setUp(self, mock_method):
        self.mock_open_file = mock_method
        self.data_path = sentinel.data_path
        self.coincidences_group = sentinel.coincidences_group

        self.cq = coincidence_queries.CoincidenceQuery(self.data_path, self.coincidences_group)

    def test_init_opens_file_and_gets_nodes(self):
        self.mock_open_file.assert_called_once_with(self.data_path, 'r')
        expected = [call(self.coincidences_group, 'coincidences'),
                    call(self.coincidences_group, 's_index'),
                    call(self.coincidences_group, 'c_index')]
        self.assertEqual(self.mock_open_file.return_value.get_node.call_args_list, expected)

        self.assertEqual(self.cq.coincidences, self.mock_open_file.return_value.get_node.return_value)
        self.assertEqual(self.cq.s_index, self.mock_open_file.return_value.get_node.return_value)
        self.assertEqual(self.cq.c_index, self.mock_open_file.return_value.get_node.return_value)

    def test_all_coincidences(self):
        result = self.cq.all_coincidences()
        self.cq.coincidences.read.assert_called_once_with()
        self.assertEqual(result, self.cq.coincidences.read.return_value)

    @patch.object(coincidence_queries.CoincidenceQuery, 'perform_query')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_s_columns')
    def test_any(self, mock_columns, mock_query):
        mock_columns.return_value = ['s501', 's502']
        self.cq.any(sentinel.stations)
        mock_columns.assert_called_once_with(sentinel.stations)
        mock_query.assert_called_once_with('s501 | s502')

    @patch.object(coincidence_queries.CoincidenceQuery, 'perform_query')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_s_columns')
    def test_all(self, mock_columns, mock_query):
        mock_columns.return_value = ['s501', 's502']
        self.cq.all(sentinel.stations)
        mock_columns.assert_called_once_with(sentinel.stations)
        mock_query.assert_called_once_with('s501 & s502')

    @patch.object(coincidence_queries.CoincidenceQuery, 'perform_query')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_s_columns')
    def test_at_least(self, mock_columns, mock_query):
        mock_columns.return_value = ['s501', 's502', 's503']
        n = 2
        self.cq.at_least(sentinel.stations, n)
        mock_columns.assert_called_once_with(sentinel.stations)
        mock_query.assert_called_once_with('(s501 & s502) | (s501 & s503) | (s502 & s503)')

    def test_perform_query(self):
        result = self.cq.perform_query(sentinel.query)
        self.cq.coincidences.read_where.assert_called_once_with(sentinel.query)
        self.assertEqual(result, self.cq.coincidences.read_where.return_value)

    def test_get_s_columns(self):
        result = self.cq._get_s_columns([501])
        self.assertEqual(result, ['s501'])


if __name__ == '__main__':
    unittest.main()
