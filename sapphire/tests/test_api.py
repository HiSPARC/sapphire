import unittest

from sapphire import api


STATION = 501


class NetworkTests(unittest.TestCase):
    def setUp(self):
        self.network = api.Network()

    def test_keys(self):
        keys = ['name', 'number']
        self.assertEqual(self.network.all_countries[0].keys(), keys)
        self.assertEqual(self.network.all_clusters[0].keys(), keys)
        self.assertEqual(self.network.all_subclusters[0].keys(), keys)
        self.assertEqual(self.network.all_stations[0].keys(), keys)
        nested_network = self.network.nested_network()
        self.assertEqual(nested_network[0].keys(), ['clusters', 'name', 'number'])
        self.assertEqual(nested_network[0]['clusters'][0].keys(), ['subclusters', 'name', 'number'])
        self.assertEqual(nested_network[0]['clusters'][0]['subclusters'][0].keys(), ['stations', 'name', 'number'])

    def test_all_is_no_filter(self):
        self.assertEqual(self.network.all_clusters, self.network.clusters())
        self.assertEqual(self.network.all_subclusters, self.network.subclusters())
        self.assertEqual(self.network.all_stations, self.network.stations())

    def test_bad_country(self):
        bad_country = 1
        self.assertRaises(Exception, self.network.clusters, country=bad_country)
        self.assertRaises(Exception, self.network.subclusters, country=bad_country)
        self.assertRaises(Exception, self.network.stations, country=bad_country)

    def test_bad_cluster(self):
        bad_cluster = 1
        self.assertRaises(Exception, self.network.subclusters, cluster=bad_cluster)
        self.assertRaises(Exception, self.network.stations, cluster=bad_cluster)

    def test_bad_subcluster(self):
        bad_subcluster = 1
        self.assertRaises(Exception, self.network.stations, subcluster=bad_subcluster)


class StationTests(unittest.TestCase):
    def setUp(self):
        self.station = api.Station(STATION)

    def test_id_numbers(self):
        self.assertEqual(self.station.station, STATION)

    def test_grouping(self):
        self.assertEqual(self.station.country, 'Netherlands')
        self.assertEqual(self.station.cluster, 'Amsterdam')
        self.assertEqual(self.station.subcluster, 'Science Park')

    def test_num_events(self):
        self.assertIsInstance(self.station.n_events(2003), int)
        self.assertEqual(self.station.n_events(2003), 0)
        self.assertEqual(self.station.n_events(2013, 8, 1), 63735)

    def test_detectors(self):
        self.assertEqual(self.station.n_detectors, 4)

    def test_event_trace(self):
        self.assertEqual(self.station.event_trace(1378771205, 571920029)[3][9], 268)


if __name__ == '__main__':
    unittest.main()
