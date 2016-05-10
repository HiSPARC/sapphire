import unittest

from mock import sentinel, patch, call

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
                    call(self.coincidences_group, 'c_index'),
                    call(self.coincidences_group, 's_index'),
                    call(self.coincidences_group, 'reconstructions')]
        call_list = self.mock_open_file.return_value.get_node.call_args_list
        self.assertEqual(call_list, expected)

        a_node = self.mock_open_file.return_value.get_node.return_value
        self.assertEqual(self.cq.coincidences, a_node)
        self.assertEqual(self.cq.s_index, a_node)
        self.assertEqual(self.cq.c_index, a_node)

    def test_finish(self):
        self.cq.finish()
        self.mock_open_file.return_value.close.assert_called_once_with()

    def test_all_coincidences(self):
        result = self.cq.all_coincidences()
        self.cq.coincidences.read.assert_called_once_with()
        self.assertEqual(result, self.cq.coincidences.read.return_value)

        result = self.cq.all_coincidences(iterator=True)
        self.cq.coincidences.iterrows.assert_called_once_with()
        self.assertEqual(result, self.cq.coincidences.iterrows.return_value)

    @patch.object(coincidence_queries.CoincidenceQuery, 'perform_query')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_allowed_s_columns')
    def test_any(self, mock_columns, mock_query):
        mock_columns.return_value = ['s501', 's502']
        self.cq.any(sentinel.stations)
        mock_columns.assert_called_once_with(sentinel.stations)
        mock_query.assert_called_once_with('(s501 | s502)', False)

    @patch.object(coincidence_queries.CoincidenceQuery, '_add_timestamp_filter')
    @patch.object(coincidence_queries.CoincidenceQuery, 'perform_query')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_allowed_s_columns')
    def test_all(self, mock_columns, mock_query, mock_ts_filter):
        mock_columns.return_value = ['s501', 's502']
        mock_ts_filter.return_value = sentinel.query
        self.cq.all([sentinel.station1, sentinel.station2])
        mock_columns.assert_called_once_with([sentinel.station1, sentinel.station2])
        mock_ts_filter.assert_called_once_with('(s501 & s502)', None, None)
        mock_query.assert_called_once_with(sentinel.query, False)
        self.assertEqual(self.cq.all([sentinel.station1, sentinel.station2, sentinel.station3]), [])

    @patch.object(coincidence_queries.CoincidenceQuery, 'perform_query')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_allowed_s_columns')
    def test_at_least(self, mock_columns, mock_query):
        mock_columns.return_value = ['s501', 's502', 's503']
        n = 2
        self.cq.at_least(sentinel.stations, n)
        mock_columns.assert_called_once_with(sentinel.stations)
        mock_query.assert_called_once_with('((s501 & s502) | (s501 & s503) | '
                                           '(s502 & s503))', False)

    @patch.object(coincidence_queries.CoincidenceQuery, 'perform_query')
    def test_timerange(self, mock_query):
        mock_query.return_value = sentinel.coincidences
        result = self.cq.timerange(1, 2)
        mock_query.assert_called_once_with('(1 <= timestamp) & (timestamp < 2)', False)
        self.assertEqual(result, sentinel.coincidences)

    def test__add_timestamp_filter(self):
        result = self.cq._add_timestamp_filter(sentinel.query)
        self.assertEqual(result, sentinel.query)
        result = self.cq._add_timestamp_filter('[query]', 1, 2)
        self.assertEqual(result, '[query] & (1 <= timestamp) & (timestamp < 2)')
        result = self.cq._add_timestamp_filter('[query]', 1)
        self.assertEqual(result, '[query] & (1 <= timestamp)')
        result = self.cq._add_timestamp_filter('[query]', stop=2)
        self.assertEqual(result, '[query] & (timestamp < 2)')

    def test_perform_query(self):
        result = self.cq.perform_query(sentinel.query)
        self.cq.coincidences.read_where.assert_called_once_with(sentinel.query)
        self.assertEqual(result, self.cq.coincidences.read_where.return_value)

    @patch.object(coincidence_queries.CoincidenceQuery, '_get_s_columns')
    def test__get_allowed_s_columns(self, mock_columns):
        self.cq.coincidences.colnames = [sentinel.scolumn1]
        mock_columns.return_value = [sentinel.scolumn1, sentinel.scolumn2]
        result = self.cq._get_allowed_s_columns([sentinel.station1, sentinel.station2])
        mock_columns.assert_called_once_with([sentinel.station1, sentinel.station2])
        self.assertEqual(result, set([sentinel.scolumn1]))

    def test__get_s_columns(self):
        result = self.cq._get_s_columns([501])
        self.assertEqual(result, ['s501'])

    @unittest.skip('WIP')
    def test__get_events(self):
        pass

    @patch.object(coincidence_queries.CoincidenceQuery, '_get_events')
    def test_all_events(self, mock_get_events):
        mock_get_events.return_value = [sentinel.event1, sentinel.event2]
        coincidences = [sentinel.coincidence]
        result = list(self.cq.all_events(coincidences))
        mock_get_events.assert_called_once_with(sentinel.coincidence)
        self.assertEqual(result, [[sentinel.event1, sentinel.event2]])

    def test_minimum_events_for_coincidence(self):
        coincidences_events = [[1], [2, 2], [3, 3, 3], [4, 4, 4, 4]]
        filtered = self.cq.minimum_events_for_coincidence(coincidences_events, 0)
        self.assertEqual(list(filtered), coincidences_events)
        filtered = self.cq.minimum_events_for_coincidence(coincidences_events)
        self.assertEqual(list(filtered), coincidences_events[1:])
        filtered = self.cq.minimum_events_for_coincidence(coincidences_events, 5)
        self.assertEqual(list(filtered), [])

    @patch.object(coincidence_queries.CoincidenceQuery, 'minimum_events_for_coincidence')
    @patch.object(coincidence_queries.CoincidenceQuery, '_events_from_stations')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_events')
    def test_events_from_stations(self, mock_get_events, mock_events_from, mock_minimum):
        mock_get_events.return_value = sentinel.events
        mock_events_from.return_value = sentinel.coincidence_events
        coincidences = [sentinel.coincidence]
        self.cq.events_from_stations(coincidences, sentinel.stations)
        # mock_get_events.assert_called_once_with(sentinel.coincidence)
        # mock_events_from.assert_called_once_with(sentinel.events, sentinel.stations)
        # mock_minimum.assert_called_once_with([sentinel.coincidence_events])

    @patch.object(coincidence_queries.CoincidenceQuery, 'minimum_events_for_coincidence')
    @patch.object(coincidence_queries.CoincidenceQuery, '_events_from_stations')
    @patch.object(coincidence_queries.CoincidenceQuery, '_get_reconstructions')
    def test_reconstructions_from_stations(self, mock_get_reconstructions, mock_events_from, mock_minimum):
        mock_get_reconstructions.return_value = sentinel.reconstructions
        mock_events_from.return_value = sentinel.coincidence_reconstructions
        coincidences = [sentinel.coincidence]
        self.cq.reconstructions_from_stations(coincidences, sentinel.stations)
        # mock_get_events.assert_called_once_with(sentinel.coincidence)
        # mock_events_from.assert_called_once_with(sentinel.events, sentinel.stations)
        # mock_minimum.assert_called_once_with([sentinel.coincidence_events])

    def test__events_from_stations(self):
        events = ([sentinel.station1, sentinel.event1],
                  [sentinel.station2, sentinel.event2])
        stations = [sentinel.station2]
        result = self.cq._events_from_stations(events, stations)
        self.assertEqual(result, [[sentinel.station2, sentinel.event2]])

    @patch.object(coincidence_queries.CoincidenceQuery, 'events_from_stations')
    @patch.object(coincidence_queries.api, 'Network')
    def test_events_in_subcluster(self, mock_network, mock_events_from):
        mock_network.return_value.station_numbers.return_value = sentinel.numbers
        mock_events_from.return_value = sentinel.coincidence_events
        result = self.cq.events_in_subcluster(sentinel.coincidences, sentinel.subcluster)
        mock_network.assert_called_once_with()
        mock_network.return_value.station_numbers.assert_called_once_with(subcluster=sentinel.subcluster)
        mock_events_from.assert_called_once_with(sentinel.coincidences, sentinel.numbers, 2)
        self.assertEqual(result, sentinel.coincidence_events)

    @patch.object(coincidence_queries.CoincidenceQuery, 'events_from_stations')
    @patch.object(coincidence_queries.api, 'Network')
    def test_events_in_cluster(self, mock_network, mock_events_from):
        mock_network.return_value.station_numbers.return_value = sentinel.numbers
        mock_events_from.return_value = sentinel.coincidence_events
        result = self.cq.events_in_cluster(sentinel.coincidences, sentinel.cluster)
        mock_network.assert_called_once_with()
        mock_network.return_value.station_numbers.assert_called_once_with(cluster=sentinel.cluster)
        mock_events_from.assert_called_once_with(sentinel.coincidences, sentinel.numbers, 2)
        self.assertEqual(result, sentinel.coincidence_events)


if __name__ == '__main__':
    unittest.main()
