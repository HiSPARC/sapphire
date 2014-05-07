from mock import sentinel, Mock, patch, call
import unittest

from sapphire.analysis import coincidence_queries


class BaseCoincidenceQueryTest(unittest.TestCase):

    @patch.object(coincidence_queries.tables, 'open_file')
    def setUp(self, mock_method):
        self.mock_open_file = mock_method
        self.data_path = sentinel.data_path
        self.coincidences_group = sentinel.coincidences_group

        self.cq = coincidence_queries.CoincidenceQuery(
            self.data_path, self.coincidences_group)

    def test_init_opens_file_and_gets_nodes(self):
        self.mock_open_file.assert_called_once_with(self.data_path, 'r')
        expected = [call(self.coincidences_group, 'coincidences'),
                    call(self.coincidences_group, 's_index'),
                    call(self.coincidences_group, 'c_index')]
        call_list = self.mock_open_file.return_value.get_node.call_args_list
        self.assertEqual(call_list, expected)

        a_node = self.mock_open_file.return_value.get_node.return_value
        self.assertEqual(self.cq.coincidences, a_node)
        self.assertEqual(self.cq.s_index, a_node)
        self.assertEqual(self.cq.c_index, a_node)

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
        mock_query.assert_called_once_with('(s501 & s502) | (s501 & s503) | '
                                           '(s502 & s503)')

    def test_perform_query(self):
        result = self.cq.perform_query(sentinel.query)
        self.cq.coincidences.read_where.assert_called_once_with(sentinel.query)
        self.assertEqual(result, self.cq.coincidences.read_where.return_value)

    def test_get_s_columns(self):
        result = self.cq._get_s_columns([501])
        self.assertEqual(result, ['s501'])

    def test_minimum_events_for_coincidence(self):
        coincidences_events = [[1], [2, 2], [3, 3, 3], [4, 4, 4, 4]]
        filtered = self.cq.minimum_events_for_coincidence(coincidences_events, 0)
        self.assertEqual(filtered, coincidences_events)
        filtered = self.cq.minimum_events_for_coincidence(coincidences_events)
        self.assertEqual(filtered, coincidences_events[1:])
        filtered = self.cq.minimum_events_for_coincidence(coincidences_events, 5)
        self.assertEqual(filtered, [])

    @patch.object(coincidence_queries.CoincidenceQuery, 'minimum_events_for_coincidence')
    @patch.object(coincidence_queries.CoincidenceQuery, '_events_from_stations')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_events')
    def test_events_from_stations(self, mock_get_events, mock_events_from, mock_minimum):
        mock_get_events.return_value = sentinel.events
        mock_events_from.return_value = sentinel.coincidence_events
        coincidences = [sentinel.coincidence]
        result = self.cq.events_from_stations(coincidences, sentinel.stations)
        mock_get_events.assert_called_once_with(sentinel.coincidence)
        mock_events_from.assert_called_once_with(sentinel.events, sentinel.stations)
        mock_minimum.assert_called_once_with([sentinel.coincidence_events])


if __name__ == '__main__':
    unittest.main()
