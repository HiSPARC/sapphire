import unittest
from datetime import date
from urllib2 import HTTPError, URLError
import warnings

from mock import patch, sentinel

from sapphire import api


STATION = 501


class APITests(unittest.TestCase):
    def setUp(self):
        self.api = api.API()

    @patch.object(api, 'urlopen')
    def test_no_check_connection(self, mock_urlopen):
        self.assertTrue(self.api.check_connection())
        mock_urlopen.return_value.read.side_effect = URLError('no interwebs!')
        self.assertFalse(self.api.check_connection())

    @patch.object(api, 'urlopen')
    def test__retrieve_url(self, mock_urlopen):
        mock_urlopen.return_value.read.side_effect = HTTPError(None, None, None, None, None)
        self.assertRaises(Exception, self.api._retrieve_url, '')
        mock_urlopen.return_value.read.side_effect = URLError('no interwebs!')
        self.assertRaises(Exception, self.api._retrieve_url, '')

    @patch.object(api, 'urlopen')
    def test__get_csv(self, mock_urlopen):
        mock_urlopen.return_value.read.return_value = '1297956608\t52.3414237\t4.8807081\t43.32'
        self.assertEqual(self.api._get_csv('gps/2/', allow_stale=False).tolist(),
                         [(1297956608, 52.3414237, 4.8807081, 43.32)])
        mock_urlopen.return_value.read.side_effect = URLError('no interwebs!')
        self.assertRaises(Exception, self.api._get_csv, 'gps/2/', allow_stale=False)
        self.assertRaises(Exception, self.api._get_csv, 'gps/0/')
        with warnings.catch_warnings(record=True) as warned:
            self.assertEqual(self.api._get_csv('gps/2/', allow_stale=True).tolist()[0],
                             (1297956608, 52.3414237, 4.8807081, 43.32))
        self.assertEqual(len(warned), 1)


@unittest.skipUnless(api.API.check_connection(), "Internet connection required")
class NetworkTests(unittest.TestCase):
    def setUp(self):
        self.network = api.Network()
        self.keys = ['name', 'number']

    @patch.object(api.Network, 'countries')
    @patch.object(api.Network, 'clusters')
    @patch.object(api.Network, 'subclusters')
    @patch.object(api.Network, 'stations')
    def test_nested_network(self, mock_stations, mock_subcluster,
                            mock_clusters, mock_countries):
        mock_countries.return_value = [{'name': sentinel.country_name,
                                        'number': sentinel.country_number}]
        mock_clusters.return_value = [{'name': sentinel.cluster_name,
                                       'number': sentinel.cluster_number}]
        mock_subcluster.return_value = [{'name': sentinel.subcluster_name,
                                         'number': sentinel.subcluster_number}]
        mock_stations.return_value = [{'name': sentinel.station_name,
                                       'number': sentinel.station_number}]
        nested_network = self.network.nested_network()
        self.assertEqual(nested_network,
                         [{'clusters': [
                           {'subclusters': [
                            {'stations': [
                             {'name': sentinel.station_name,
                              'number': sentinel.station_number}],
                             'name': sentinel.subcluster_name,
                             'number': sentinel.subcluster_number}],
                            'name': sentinel.cluster_name,
                            'number': sentinel.cluster_number}],
                           'name': sentinel.country_name,
                           'number': sentinel.country_number}])

    def test_lazy_countries(self):
        self.laziness_of_method('countries')

    def test_countries(self):
        self.network.countries()
        self.assertEqual(self.network._all_countries, self.network.countries())
        self.assertEqual(self.network.countries()[0].keys(), self.keys)

    def test_lazy_clusters(self):
        self.laziness_of_method('clusters')

    def test_clusters(self):
        self.network.clusters()
        self.assertEqual(self.network._all_clusters, self.network.clusters())
        self.assertEqual(self.network.clusters()[0].keys(), self.keys)
        self.assertEqual(self.network.clusters(country=70000)[0]['number'], 70000)

    def test_bad_clusters(self):
        bad_number = 1
        self.assertRaises(Exception, self.network.clusters, country=bad_number)

    def test_lazy_subclusters(self):
        self.laziness_of_method('subclusters')

    def test_subcluster(self):
        self.network.subclusters()
        self.assertEqual(self.network._all_subclusters, self.network.subclusters())
        self.assertEqual(self.network.subclusters()[0].keys(), self.keys)
        self.assertEqual(self.network.subclusters(country=70000)[0]['number'], 70000)
        self.assertEqual(self.network.subclusters(cluster=70000)[0]['number'], 70000)

    def test_bad_subcluster(self):
        bad_number = 1
        self.assertRaises(Exception, self.network.subclusters, country=bad_number)
        self.assertRaises(Exception, self.network.subclusters, cluster=bad_number)

    def test_country_numbers(self):
        self.network._all_countries = [{'number': sentinel.number1},
                                       {'number': sentinel.number2}]
        self.assertEqual(self.network.country_numbers(),
                         [sentinel.number1, sentinel.number2])

    @patch.object(api.Network, 'clusters')
    @patch.object(api.Network, 'validate_numbers')
    def test_cluster_numbers(self, mock_validate, mock_clusters):
        mock_clusters.return_value = [{'number': sentinel.number1},
                                      {'number': sentinel.number2}]
        self.assertEqual(self.network.cluster_numbers(sentinel.country),
                         [sentinel.number1, sentinel.number2])
        mock_clusters.assert_called_once_with(country=sentinel.country)

    @patch.object(api.Network, 'subclusters')
    @patch.object(api.Network, 'validate_numbers')
    def test_subcluster_numbers(self, mock_validate, mock_subclusters):
        mock_subclusters.return_value = [{'number': sentinel.number1},
                                         {'number': sentinel.number2}]
        self.assertEqual(self.network.subcluster_numbers(sentinel.country,
                                                         sentinel.cluster),
                         [sentinel.number1, sentinel.number2])
        mock_subclusters.assert_called_once_with(country=sentinel.country,
                                                 cluster=sentinel.cluster)

    @patch.object(api.Network, 'stations')
    @patch.object(api.Network, 'validate_numbers')
    def test_station_numbers(self, mock_validate, mock_stations):
        mock_stations.return_value = [{'number': sentinel.number1},
                                      {'number': sentinel.number2}]
        self.assertEqual(self.network.station_numbers(sentinel.country,
                                                      sentinel.cluster,
                                                      sentinel.subcluster),
                         [sentinel.number1, sentinel.number2])
        mock_stations.assert_called_once_with(country=sentinel.country,
                                              cluster=sentinel.cluster,
                                              subcluster=sentinel.subcluster)

    @patch.object(api.Network, '_get_json')
    def test_no_stale_station_numbers(self, mock_get_json):
        mock_get_json.side_effect = Exception('no interwebs!')
        self.assertRaises(Exception, self.network.station_numbers, allow_stale=False)

    @patch.object(api.Network, '_get_json')
    def test_stale_station_numbers(self, mock_get_json):
        mock_get_json.side_effect = Exception('no interwebs!')
        with warnings.catch_warnings(record=True) as warned:
            station_numbers = self.network.station_numbers(allow_stale=True)
            self.assertEqual(min(station_numbers), 2)
            station_numbers = self.network.station_numbers(country=20000, allow_stale=True)
            self.assertEqual(min(station_numbers), 20001)
            station_numbers = self.network.station_numbers(cluster=1000, allow_stale=True)
            self.assertEqual(min(station_numbers), 1001)
            station_numbers = self.network.station_numbers(subcluster=500, allow_stale=True)
            self.assertEqual(min(station_numbers), 501)
        self.assertEqual(len(warned), 1)

    def test_invalid_query_for_station_numbers(self):
        bad_number = 1
        for allow in (True, False):
            self.assertRaises(Exception, self.network.station_numbers,
                              country=bad_number, allow_stale=allow)
            self.assertRaises(Exception, self.network.station_numbers,
                              cluster=bad_number, allow_stale=allow)
            self.assertRaises(Exception, self.network.station_numbers,
                              subcluster=bad_number, allow_stale=allow)

    def test_lazy_stations(self):
        self.laziness_of_method('stations')

    def test_stations(self):
        self.network.stations()
        self.assertEqual(self.network._all_stations, self.network.stations())
        self.assertEqual(self.network.stations()[0].keys(), self.keys)
        self.assertEqual(self.network.stations(country=70000)[0]['number'], 70001)
        self.assertEqual(self.network.stations(cluster=70000)[0]['number'], 70001)
        self.assertEqual(self.network.stations(subcluster=70000)[0]['number'], 70001)

    def test_bad_stations(self):
        bad_number = 1
        self.assertRaises(Exception, self.network.stations, country=bad_number)
        self.assertRaises(Exception, self.network.stations, cluster=bad_number)
        self.assertRaises(Exception, self.network.stations, subcluster=bad_number)

    def test_stations_with_data(self):
        stations_with_data = self.network.stations_with_data(2004, 1, 9)
        self.assertEqual(stations_with_data[0].keys(), self.keys)
        self.assertEqual(stations_with_data[0]['number'], 2)
        self.assertEqual(len(self.network.stations_with_data(2004, 1, 1)), 0)
        self.assertRaises(Exception, self.network.stations_with_data, year=2004, day=1)
        self.assertRaises(Exception, self.network.stations_with_data, month=1, day=1)
        self.assertRaises(Exception, self.network.stations_with_data, month=1)
        self.assertRaises(Exception, self.network.stations_with_data, day=1)

    def test_stations_with_weather(self):
        stations_with_weather = self.network.stations_with_weather(2004, 10, 9)
        self.assertEqual(stations_with_weather[0].keys(), self.keys)
        self.assertEqual(stations_with_weather[0]['number'], 3)
        self.assertEqual(len(self.network.stations_with_weather(2004, 1, 1)), 0)
        self.assertRaises(Exception, self.network.stations_with_weather, year=2004, day=1)
        self.assertRaises(Exception, self.network.stations_with_weather, month=1, day=1)
        self.assertRaises(Exception, self.network.stations_with_weather, month=1)
        self.assertRaises(Exception, self.network.stations_with_weather, day=1)

    def test_coincidence_time(self):
        names = ('hour', 'counts')
        data = self.network.coincidence_time(2013, 1, 1)
        self.assertEqual(data.dtype.names, names)
        self.assertTrue((data['hour'] == range(24)).all())
        self.assertEqual(data['counts'][0], 424)

    def test_coincidence_number(self):
        names = ('n', 'counts')
        data = self.network.coincidence_number(2013, 1, 1)
        self.assertEqual(data.dtype.names, names)
        self.assertTrue((data['n'] == range(2, 100)).all())
        self.assertEqual(data['counts'][0], 8763)

    def laziness_of_method(self, method):
        with patch.object(api.API, '_get_json') as mock_get_json:
            self.assertFalse(mock_get_json.called)
            data = self.network.__getattribute__(method)()
            self.assertTrue(mock_get_json.called)
            self.assertEqual(mock_get_json.call_count, 1)
            data2 = self.network.__getattribute__(method)()
            self.assertEqual(mock_get_json.call_count, 1)
            self.assertEqual(data, data2)


@unittest.skipUnless(api.API.check_connection(), "Internet connection required")
class StationTests(unittest.TestCase):
    def setUp(self):
        self.station = api.Station(STATION)

    @patch.object(api.Station, '_get_json')
    def test_no_stale_station(self, mock_get_json):
        mock_get_json.side_effect = Exception('no interwebs!')
        self.assertRaises(Exception, api.Station, 501, allow_stale=False)

    @patch.object(api.Station, '_get_json')
    def test_stale_station(self, mock_get_json):
        mock_get_json.side_effect = Exception('no interwebs!')
        with warnings.catch_warnings(record=True) as warned:
            api.Station(501, allow_stale=True)
        self.assertEqual(len(warned), 1)

    def test_id_numbers(self):
        self.assertEqual(self.station.station, STATION)

    def test_properties(self):
        self.assertEqual(self.station.country(), 'Netherlands')
        self.assertEqual(self.station.cluster(), 'Amsterdam')
        self.assertEqual(self.station.subcluster(), 'Science Park')
        self.assertEqual(self.station.n_detectors(), 4)

    def test_detectors(self):
        keys = ['mpv', 'alpha', 'beta', 'radius', 'height']
        self.assertEqual(len(self.station.detectors()), self.station.n_detectors())
        self.assertEqual(self.station.detectors()[0].keys(), keys)
        self.assertEqual(self.station.detectors(date(2011, 1, 1))[0].keys(), keys)
        self.assertEqual(self.station.detectors(date(2011, 1, 1))[0]['alpha'], 225)

    def test_location(self):
        keys = ['latitude', 'altitude', 'longitude']
        self.assertEqual(self.station.location().keys(), keys)
        self.assertEqual(self.station.location(date(2002, 1, 1))['latitude'],
                         52.3559179545407)

    def test_config(self):
        self.assertEqual(self.station.config()['detnum'], 501)
        self.assertEqual(self.station.config(date(2011, 1, 1))['mas_ch1_current'],
                         7.54901960784279)

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

    @patch.object(api, 'urlopen')
    def test_event_trace(self, mock_urlopen):
        trace = '[%s]' % ', '.join(str(v) for v in range(0, 11))
        mock_urlopen.return_value.read.return_value = '[%s]' % ', '.join(4 * [trace])
        self.assertEqual(self.station.event_trace(1378771205, 571920029)[3][9], 9)
        trace = '[%s]' % ', '.join(str(v) for v in range(200, 211))
        mock_urlopen.return_value.read.return_value = '[%s]' % ', '.join(4 * [trace])
        self.assertEqual(self.station.event_trace(1378771205, 571920029, raw=True)[3][9], 209)

    def test_event_time(self):
        names = ('hour', 'counts')
        data = self.station.event_time(2013, 1, 1)
        self.assertEqual(data.dtype.names, names)
        self.assertTrue((data['hour'] == range(24)).all())

    def test_pulse_height(self):
        names = ('pulseheight', 'ph1', 'ph2', 'ph3', 'ph4')
        data = self.station.pulse_height(2013, 1, 1)
        self.assertEqual(data.dtype.names, names)
        self.assertTrue((data['pulseheight'] == range(0, 2500, 10)).all())

    def test_pulse_integral(self):
        names = ('pulseintegral', 'pi1', 'pi2', 'pi3', 'pi4')
        data = self.station.pulse_integral(2013, 1, 1)
        self.assertEqual(data.dtype.names, names)
        self.assertTrue((data['pulseintegral'] == range(0, 62500, 250)).all())

    def test_barometer(self):
        names = ('timestamp', 'air_pressure')
        data = self.station.barometer(2013, 1, 1)
        self.assertEqual(data.dtype.names, names)

    def test_temperature(self):
        names = ('timestamp', 'temperature')
        data = self.station.temperature(2013, 1, 1)
        self.assertEqual(data.dtype.names, names)

    def test_voltages(self):
        names = ('timestamp', 'voltage1', 'voltage2', 'voltage3', 'voltage4')
        data = self.station.voltages
        self.assertEqual(data.dtype.names, names)

    def test_laziness_voltages(self):
        self.laziness_of_attribute('voltages')

    def test_voltage(self):
        data = self.station.voltage(1378771200)  # 2013-9-10
        self.assertEqual(data, [954, 860, 714, 752])

        data = self.station.voltage(0)  # 1970-1-1
        data2 = self.station.voltages[0]
        self.assertEqual(data, [data2['voltage1'], data2['voltage2'],
                                data2['voltage3'], data2['voltage4']])
        data = self.station.voltage(2208988800)  # 2040-1-1
        data2 = self.station.voltages[-1]
        self.assertEqual(data, [data2['voltage1'], data2['voltage2'],
                                data2['voltage3'], data2['voltage4']])

    def test_laziness_currents(self):
        self.laziness_of_attribute('currents')

    def test_currents(self):
        names = ('timestamp', 'current1', 'current2', 'current3', 'current4')
        data = self.station.currents
        self.assertEqual(data.dtype.names, names)

    def test_current(self):
        data = self.station.current(1378771200)  # 2013-9-10
        self.assertEqual(data, [7.84, 7.94, 10.49, 10.88])

    def test_laziness_gps_locations(self):
        self.laziness_of_attribute('gps_locations')

    def test_gps_locations(self):
        names = ('timestamp', 'latitude', 'longitude', 'altitude')
        data = self.station.gps_locations
        self.assertEqual(data.dtype.names, names)

    def test_gps_location(self):
        keys = ['latitude', 'longitude', 'altitude']
        data = self.station.gps_location(1378771200)  # 2013-9-10
        self.assertItemsEqual(data.keys(), keys)
        self.assertItemsEqual(data.values(), [52.3559286, 4.9511443, 54.97])

    def test_laziness_station_layouts(self):
        self.laziness_of_attribute('station_layouts')

    def test_station_layouts(self):
        names = ('timestamp',
                 'radius1', 'alpha1', 'height1', 'beta1',
                 'radius2', 'alpha2', 'height2', 'beta2',
                 'radius3', 'alpha3', 'height3', 'beta3',
                 'radius4', 'alpha4', 'height4', 'beta4')
        data = self.station.station_layouts
        self.assertEqual(data.dtype.names, names)

    def test_station_layout(self):
        data = self.station.station_layout(0)
        self.assertEqual(len(data), 4)
        self.assertEqual(len(data[0]), 4)

    def test_laziness_detector_timing_offsets(self):
        self.laziness_of_attribute('detector_timing_offsets')

    @patch.object(api, 'urlopen')
    def test_detector_timing_offsets(self, mock_urlopen):
        mock_urlopen.return_value.read.return_value = '1234567980\t0.0\t2.5\t-2.5\t0.25\n' * 4
        names = ('timestamp', 'offset1', 'offset2', 'offset3', 'offset4')
        data = self.station.detector_timing_offsets
        self.assertEqual(data.dtype.names, names)
        self.assertEqual(len(data), 4)
        self.assertEqual(len(data[0]), 5)

    @patch.object(api, 'urlopen')
    def test_detector_timing_offset(self, mock_urlopen):
        mock_urlopen.return_value.read.return_value = '1234567980\t0.0\t2.5\t-2.5\t0.25\n' * 4
        offsets = self.station.detector_timing_offset(0)
        self.assertEqual(len(offsets), 4)

    def laziness_of_attribute(self, attribute):
        with patch.object(api.API, '_get_csv') as mock_get_csv:
            self.assertFalse(mock_get_csv.called)
            data = self.station.__getattribute__(attribute)
            self.assertTrue(mock_get_csv.called)
            self.assertEqual(mock_get_csv.call_count, 1)
            data2 = self.station.__getattribute__(attribute)
            self.assertEqual(mock_get_csv.call_count, 1)
            self.assertEqual(data, data2)


if __name__ == '__main__':
    unittest.main()
