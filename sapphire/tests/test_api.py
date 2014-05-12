import unittest
from datetime import date

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
        self.assertEqual(self.network.stations_with_data(2004, 1, 9)[0].keys(), keys)
        self.assertEqual(self.network.stations_with_weather(2004, 10, 9)[0].keys(), keys)
        nested_network = self.network.nested_network()
        self.assertEqual(nested_network[0].keys(), ['clusters', 'name', 'number'])
        self.assertEqual(nested_network[0]['clusters'][0].keys(), ['subclusters', 'name', 'number'])
        self.assertEqual(nested_network[0]['clusters'][0]['subclusters'][0].keys(), ['stations', 'name', 'number'])

    def test_countries(self):
        self.assertEqual(self.network.all_countries, self.network.countries())

    def test_clusters(self):
        self.assertEqual(self.network.all_clusters, self.network.clusters())
        self.assertEqual(self.network.clusters(country=70000)[0]['number'], 70000)
        bad_number = 1
        self.assertRaises(Exception, self.network.clusters, country=bad_number)

    def test_subcluster(self):
        self.assertEqual(self.network.all_subclusters, self.network.subclusters())
        self.assertEqual(self.network.subclusters(country=70000)[0]['number'], 70000)
        self.assertEqual(self.network.subclusters(cluster=70000)[0]['number'], 70000)
        bad_number = 1
        self.assertRaises(Exception, self.network.subclusters, country=bad_number)
        self.assertRaises(Exception, self.network.subclusters, cluster=bad_number)

    def test_bad_stations(self):
        self.assertEqual(self.network.all_stations, self.network.stations())
        self.assertEqual(self.network.stations(country=70000)[0]['number'], 70001)
        self.assertEqual(self.network.stations(cluster=70000)[0]['number'], 70001)
        self.assertEqual(self.network.stations(subcluster=70000)[0]['number'], 70001)
        bad_number = 1
        self.assertRaises(Exception, self.network.stations, country=bad_number)
        self.assertRaises(Exception, self.network.stations, cluster=bad_number)
        self.assertRaises(Exception, self.network.stations, subcluster=bad_number)

    def test_stations_with_data(self):
        self.assertEqual(self.network.stations_with_data(2004, 1, 9)[0]['number'], 2)
        self.assertEqual(len(self.network.stations_with_data(2004, 1, 1)), 0)
        self.assertRaises(Exception, self.network.stations_with_data, day=1)
        self.assertRaises(Exception, self.network.stations_with_data, month=1)
        self.assertRaises(Exception, self.network.stations_with_data, day=1, month=1)

    def test_stations_with_weather(self):
        self.assertEqual(self.network.stations_with_weather(2004, 10, 9)[0]['number'], 3)
        self.assertEqual(len(self.network.stations_with_weather(2004, 1, 1)), 0)
        self.assertRaises(Exception, self.network.stations_with_weather, day=1)
        self.assertRaises(Exception, self.network.stations_with_weather, month=1)
        self.assertRaises(Exception, self.network.stations_with_weather, day=1, month=1)


class StationTests(unittest.TestCase):
    def setUp(self):
        self.station = api.Station(STATION)

    def test_id_numbers(self):
        self.assertEqual(self.station.station, STATION)

    def test_properties(self):
        self.assertEqual(self.station.country, 'Netherlands')
        self.assertEqual(self.station.cluster, 'Amsterdam')
        self.assertEqual(self.station.subcluster, 'Science Park')
        self.assertEqual(self.station.n_detectors, 4)

    def test_detectors(self):
        keys = ['mpv', 'alpha', 'beta', 'radius', 'height']
        self.assertEqual(len(self.station.detectors()), self.station.n_detectors)
        self.assertEqual(self.station.detectors()[0].keys(), keys)
        self.assertEqual(self.station.detectors(date(2011, 1, 1))[0].keys(), keys)
        self.assertEqual(self.station.detectors()[0]['alpha'], 225)
        self.assertEqual(self.station.detectors(date(2011, 1, 1))[0]['alpha'], 225)

    def test_location(self):
        keys = ['latitude',  'altitude', 'longitude']
        self.assertEqual(self.station.location().keys(), keys)
        self.assertEqual(self.station.location()['latitude'], 52.355928561847)
        self.assertEqual(self.station.location(date(2002, 1, 1))['latitude'], 52.3559179545407)

    def test_config(self):
        self.assertEqual(self.station.config()['detnum'], 501)
        self.assertEqual(self.station.config(date(2011, 1, 1))['mas_ch1_current'], 7.54901960784279)

    def test_num_events(self):
        self.assertIsInstance(self.station.n_events(2003), int)
        self.assertEqual(self.station.n_events(2003), 0)
        self.assertEqual(self.station.n_events(2013, 8, 1), 63735)
        # No year
        self.assertRaises(Exception, self.station.n_events, month=1, day=1, hour=1)
        self.assertRaises(Exception, self.station.n_events, month=1, day=1)
        self.assertRaises(Exception, self.station.n_events, month=1, hour=1)
        self.assertRaises(Exception, self.station.n_events, month=1)
        self.assertRaises(Exception, self.station.n_events, day=1, hour=1)
        self.assertRaises(Exception, self.station.n_events, day=1)
        self.assertRaises(Exception, self.station.n_events, hour=1)
        # No month
        self.assertRaises(Exception, self.station.n_events, year=2011, day=1, hour=1)
        self.assertRaises(Exception, self.station.n_events, year=2011, day=1)
        self.assertRaises(Exception, self.station.n_events, year=2011, hour=1)
        # No day
        self.assertRaises(Exception, self.station.n_events, year=2011, month=1, hour=1)

    def test_has_data(self):
        self.assertEqual(self.station.has_data(), True)
        self.assertEqual(self.station.has_data(2014), True)
        self.assertEqual(self.station.has_data(2002, 1), False)
        self.assertEqual(self.station.has_data(2014, 1, 1), True)
        self.assertRaises(Exception, self.station.has_data, day=1)
        self.assertRaises(Exception, self.station.has_data, month=1)
        self.assertRaises(Exception, self.station.has_data, month=1, day=1)
        self.assertRaises(Exception, self.station.has_data, year=2011, day=1)

    def test_has_weather(self):
        self.assertEqual(self.station.has_weather(), True)
        self.assertEqual(self.station.has_weather(2014), True)
        self.assertEqual(self.station.has_weather(2002, 1), False)
        self.assertEqual(self.station.has_weather(2014, 1, 1), True)
        self.assertRaises(Exception, self.station.has_weather, day=1)
        self.assertRaises(Exception, self.station.has_weather, month=1)
        self.assertRaises(Exception, self.station.has_weather, month=1, day=1)
        self.assertRaises(Exception, self.station.has_weather, year=2011, day=1)

    def test_event_trace(self):
        self.assertEqual(self.station.event_trace(1378771205, 571920029)[3][9], 268)


if __name__ == '__main__':
    unittest.main()
