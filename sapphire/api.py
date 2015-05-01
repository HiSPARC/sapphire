""" HiSPARC api interface

    This provides easy classes and functions to access the HiSPARC
    publicdb API. This takes care of the url retrieval and conversion
    from JSON to Python dictionaries.

    Example usage:

    .. code-block:: python

        >>> from sapphire.api import Station
        >>> stations = [5, 301, 3102, 504, 7101, 8008, 13005]
        >>> clusters = [Station(station).cluster.lower()
        ...             for station in stations]
        >>> station_groups = ['/hisparc/cluster_%s/station_%d' % (c, s)
        ...                   for c, s in zip(clusters, stations)]
        >>> station_groups
        [u'/hisparc/cluster_amsterdam/station_5',
         u'/hisparc/cluster_amsterdam/station_301',
         u'/hisparc/cluster_leiden/station_3102',
         u'/hisparc/cluster_amsterdam/station_504',
         u'/hisparc/cluster_enschede/station_7101',
         u'/hisparc/cluster_eindhoven/station_8008',
         u'/hisparc/cluster_bristol/station_13005']

"""
import logging
import datetime
import json
import warnings
from os import path
from urllib2 import urlopen, HTTPError, URLError
from StringIO import StringIO
from bisect import bisect_right

from lazy import lazy
from numpy import genfromtxt, atleast_1d

logger = logging.getLogger('api')

API_BASE = 'http://data.hisparc.nl/api/'
SRC_BASE = 'http://data.hisparc.nl/show/source/'
JSON_FILE = path.join(path.dirname(__file__), 'data/hisparc_stations.json')


class API(object):

    """Base API class

    This provided the methods to retrieve data from the API. The results
    are converted from JSON data to Python objects (dict/list/etc).
    Support is also provided for the retrieval of Source CSV data, which
    is returned as NumPy arrays.

    """

    urls = {"stations": 'stations/',
            "stations_in_subcluster": 'subclusters/{subcluster_number}/',
            "subclusters": 'subclusters/',
            "subclusters_in_cluster": 'clusters/{cluster_number}/',
            "clusters": 'clusters/',
            "clusters_in_country": 'countries/{country_number}/',
            "countries": 'countries/',
            "stations_with_data": 'stations/data/{year}/{month}/{day}/',
            "stations_with_weather": 'stations/weather/{year}/{month}/{day}/',
            "station_info": 'station/{station_number}/{year}/{month}/{day}/',
            "has_data": 'station/{station_number}/data/{year}/{month}/{day}/',
            "has_weather": 'station/{station_number}/weather/{year}/{month}/'
                           '{day}/',
            "configuration": 'station/{station_number}/config/{year}/'
                             '{month}/{day}/',
            "number_of_events": 'station/{station_number}/num_events/{year}/'
                                '{month}/{day}/{hour}/',
            "event_trace": 'station/{station_number}/trace/{ext_timestamp}/',
            "pulseheight_fit": 'station/{station_number}/plate/{plate_number}/'
                               'pulseheight/fit/{year}/{month}/{day}/',
            "pulseheight_drift": 'station/{station_number}/plate/'
                                 '{plate_number}/pulseheight/drift/{year}/'
                                 '{month}/{day}/{number_of_days}/'}

    src_urls = {
        'coincidencetime': 'coincidencetime/{year}/{month}/{day}/',
        'coincidencenumber': 'coincidencenumber/{year}/{month}/{day}/',
        'eventtime': 'eventtime/{station_number}/{year}/{month}/{day}/',
        'pulseheight': 'pulseheight/{station_number}/{year}/{month}/{day}/',
        'integral': 'pulseintegral/{station_number}/{year}/{month}/{day}/',
        'barometer': 'barometer/{station_number}/{year}/{month}/{day}/',
        'temperature': 'temperature/{station_number}/{year}/{month}/{day}/',
        'voltage': 'voltage/{station_number}/',
        'current': 'current/{station_number}/',
        'gps': 'gps/{station_number}/'}

    @classmethod
    def _get_json(cls, urlpath):
        """Retrieve a JSON from the HiSPARC API

        :param urlpath: the api urlpath (after http://data.hisparc.nl/api/)
            to retrieve.
        :return: the data returned by the api as dictionary or integer.

        """
        json_data = cls._retrieve_url(urlpath)
        data = json.loads(json_data)

        return data

    @classmethod
    def _get_csv(cls, urlpath, names=None):
        """Retrieve a Source CSV from the HiSPARC Public Database

        :param urlpath: the csv urlpath to retrieve
            (after http://data.hisparc.nl/show/source/).
        :return: the data returned as array.

        """
        csv_data = cls._retrieve_url(urlpath, base=SRC_BASE)
        data = genfromtxt(StringIO(csv_data), delimiter='\t', dtype=None,
                          names=names)

        return atleast_1d(data)

    @staticmethod
    def _retrieve_url(urlpath, base=API_BASE):
        """Open a HiSPARC API URL and read the data

        :param urlpath: the api urlpath (after http://data.hisparc.nl/api/)
            to retrieve
        :return: the data returned by the api as a string

        """
        url = base + urlpath
        logging.debug('Getting: ' + url)
        try:
            result = urlopen(url).read()
        except HTTPError, e:
            raise Exception('A HTTP %d error occured for the url: %s' %
                            (e.code, url))
        except URLError:
            raise Exception('An URL error occured.')

        return result

    @staticmethod
    def check_connection():
        """Open the API man page URL to test the connection

        :return: boolean indicating the internet status

        """
        try:
            urlopen(API_BASE).read()
        except URLError:
            return False
        return True

    @staticmethod
    def get_active_index(timestamps, timestamp):
        """Get the index where the timestamp fits.

        :param timestamps: list of timestamps.
        :param timestamp: timestamp for which to find the position.
        :return: index into the timestamps list.

        """
        idx = bisect_right(timestamps, timestamp, lo=0)
        if idx == 0:
            idx = 1
        return idx - 1


class Network(API):
    """Get info about the network (countries/clusters/subclusters/stations)"""

    _all_countries = None
    _all_clusters = None
    _all_subclusters = None
    _all_stations = None

    def __init__(self):
        """Initialize network

        """
        pass

    def countries(self):
        """Get a list of countries

        :return: all countries in the region

        """
        if not self._all_countries:
            path = self.urls['countries']
            self._all_countries = self._get_json(path)
        return self._all_countries

    def country_numbers(self):
        """Same as countries but only retuns a list of country numbers"""

        countries = self.countries()
        return [country['number'] for country in countries]

    def clusters(self, country=None):
        """Get a list of clusters

        :param country: the number of the country for which to get all
            clusters.
        :return: all clusters in the region

        """
        self.validate_numbers(country)
        if country is None:
            if not self._all_clusters:
                path = self.urls['clusters']
                self._all_clusters = self._get_json(path)
            clusters = self._all_clusters
        else:
            path = (self.urls['clusters_in_country']
                    .format(country_number=country))
            clusters = self._get_json(path)
        return clusters

    def cluster_numbers(self, country=None):
        """Same as clusters but only retuns a list of cluster numbers"""

        self.validate_numbers(country)
        clusters = self.clusters(country=country)
        return [cluster['number'] for cluster in clusters]

    def subclusters(self, country=None, cluster=None):
        """Get a list of subclusters

        :param country,cluster: the number of the region for which to get
            the subclusters it contains, only one or none should
            be specified.
        :return: all subclusters in the region

        """
        self.validate_numbers(country, cluster)
        if country is None and cluster is None:
            if not self._all_subclusters:
                path = self.urls['subclusters']
                self._all_subclusters = self._get_json(path)
            subclusters = self._all_subclusters
        elif country is not None:
            subclusters = []
            path = (self.urls['clusters_in_country']
                    .format(country_number=country))
            clusters = self._get_json(path)
            for cluster in clusters:
                path = (self.urls['subclusters_in_cluster']
                        .format(cluster_number=cluster['number']))
                subclusters.extend(self._get_json(path))
        else:
            path = (self.urls['subclusters_in_cluster']
                    .format(cluster_number=cluster))
            subclusters = self._get_json(path)
        return subclusters

    def subcluster_numbers(self, country=None, cluster=None):
        """Same as subclusters but only retuns a list of subcluster numbers"""

        self.validate_numbers(country, cluster)
        subclusters = self.subclusters(country=country, cluster=cluster)
        return [subcluster['number'] for subcluster in subclusters]

    def stations(self, country=None, cluster=None, subcluster=None):
        """Get a list of stations

        :param country,cluster,subcluster: the number of the region
            for which to get all stations, only one or none should
            be specified.
        :return: all stations in the region

        """
        self.validate_numbers(country, cluster, subcluster)
        if country is None and cluster is None and subcluster is None:
            if not self._all_stations:
                path = self.urls['stations']
                self._all_stations = self._get_json(path)
            stations = self._all_stations
        elif country is not None:
            stations = []
            path = (self.urls['clusters_in_country']
                    .format(country_number=country))
            clusters = self._get_json(path)
            for cluster in clusters:
                path = (self.urls['subclusters_in_cluster']
                        .format(cluster_number=cluster['number']))
                subclusters = self._get_json(path)
                for subcluster in subclusters:
                    path = (self.urls['stations_in_subcluster']
                            .format(subcluster_number=subcluster['number']))
                    stations.extend(self._get_json(path))
        elif cluster is not None:
            stations = []
            path = (self.urls['subclusters_in_cluster']
                    .format(cluster_number=cluster))
            subclusters = self._get_json(path)
            for subcluster in subclusters:
                path = (self.urls['stations_in_subcluster']
                        .format(subcluster_number=subcluster['number']))
                stations.extend(self._get_json(path))
        else:
            path = (self.urls['stations_in_subcluster']
                    .format(subcluster_number=subcluster))
            stations = self._get_json(path)
        return stations

    def station_numbers(self, country=None, cluster=None, subcluster=None,
                        allow_stale=True):
        """Same as stations but only retuns a list of station numbers"""

        self.validate_numbers(country, cluster, subcluster)
        try:
            stations = self.stations(country=country, cluster=cluster,
                                     subcluster=subcluster)
            return [station['number'] for station in stations]
        except Exception, e:
            if not allow_stale:
                raise
            # Try getting the station info from the JSON.
            try:
                if country is None and cluster is None and subcluster is None:
                    start, end = (0, 1e9)
                elif country is not None:
                    start, end = (country, country + 10000)
                elif cluster is not None:
                    start, end = (cluster, cluster + 1000)
                else:
                    start, end = (subcluster, subcluster + 100)
                with open(JSON_FILE) as data:
                    stations = [int(s) for s in json.load(data).keys()
                                if s != '_info' and start <= int(s) < end]
                warnings.warn('Couldnt get values from the server, using '
                              'hard-coded values. Possibly outdated.',
                              UserWarning)
                return sorted(stations)
            except:
                raise e

    def nested_network(self):
        """Get a nested list of the full network"""
        countries = self.countries()
        for country in countries:
            clusters = self.clusters(country=country['number'])
            for cluster in clusters:
                subclusters = self.subclusters(cluster=cluster['number'])
                for subcluster in subclusters:
                    stations = self.stations(subcluster=subcluster['number'])
                    subcluster.update({'stations': stations})
                cluster.update({'subclusters': subclusters})
            country.update({'clusters': clusters})
        return countries

    @classmethod
    def stations_with_data(cls, year='', month='', day=''):
        """Get a list of stations with data on the specified date

        :param year,month,day: the date for which to check. It is
            possible to be less specific.
        :return: all stations with data.

        """
        if year == '' and (month != '' or day != ''):
            raise Exception('You must also specify the year')
        elif month == '' and day != '':
            raise Exception('You must also specify the month')

        path = (cls.urls['stations_with_data']
                .format(year=year, month=month, day=day).strip("/"))
        return cls._get_json(path)

    @classmethod
    def stations_with_weather(cls, year='', month='', day=''):
        """Get a list of stations with weather data on the specified date

        :param year,month,day: the date for which to check. It is
            possible to be less specific.
        :return: all stations with weather data.

        """
        if year == '' and (month != '' or day != ''):
            raise Exception('You must also specify the year')
        elif month == '' and day != '':
            raise Exception('You must also specify the month')

        path = (cls.urls['stations_with_weather']
                .format(year=year, month=month, day=day).strip("/"))
        return cls._get_json(path)

    @classmethod
    def coincidence_time(cls, year, month, day):
        """Get the coincidences per hour histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('hour', 'counts')
        path = cls.src_urls['coincidencetime'].format(year=year, month=month,
                                                      day=day)
        return cls._get_csv(path, names=columns)

    @classmethod
    def coincidence_number(cls, year, month, day):
        """Get the number of stations in coincidence histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('n', 'counts')
        path = cls.src_urls['coincidencenumber'].format(year=year, month=month,
                                                        day=day)
        return cls._get_csv(path, names=columns)

    @staticmethod
    def validate_numbers(country=None, cluster=None, subcluster=None):
        if country is not None and country % 10000:
            raise Exception('Invalid country number, '
                            'must be multiple of 10000.')
        if cluster is not None and cluster % 1000:
            raise Exception('Invalid cluster number, '
                            'must be multiple of 1000.')
        if subcluster is not None and subcluster % 100:
            raise Exception('Invalid subcluster number, '
                            'must be multiple of 100.')


class Station(API):
    """Access data about a single station"""

    def __init__(self, station, date=None, allow_stale=True):
        """Initialize station

        :param station: station number.
        :param date: date object for which to get the station information.
        :param allow_stale: set to False to require data to be fresh
                            from the server.

        """
        self.station = station
        if date is None:
            date = datetime.date.today()
        path = (self.urls['station_info']
                .format(station_number=self.station, year=date.year,
                        month=date.month, day=date.day))
        try:
            self.info = self._get_json(path)
        except Exception, e:
            if allow_stale:
                # Try getting the station info from the JSON.
                try:
                    with open(JSON_FILE) as data:
                        self.info = json.load(data)[str(station)]
                    warnings.warn('Couldnt get values from the server, using '
                                  'hard-coded values. Not all info available.',
                                  UserWarning)
                except:
                    raise e
            else:
                raise

    def country(self):
        return self.info['country']

    def cluster(self):
        return self.info['cluster']

    def subcluster(self):
        return self.info['subcluster']

    def n_detectors(self):
        """Get the number of detectors in this station"""
        return len(self.info['scintillators'])

    def detectors(self, date=None):
        """Get the locations of detectors

        :param date: date object for which to get the detector information
        :return: list with the locations of each detector

        """
        if date is None:
            return self.info['scintillators']
        else:
            station = Station(self.station, date)
            return station.detectors()

    def location(self, date=None):
        """Get gps location of the station

        :param date: date object for which to get the location
        :return: the gps coordinates for the station

        """
        if date is None:
            dict = self.info
            return {'latitude': dict['latitude'],
                    'longitude': dict['longitude'],
                    'altitude': dict['altitude']}
        else:
            dict = self.config(date=date)
            return {'latitude': dict['gps_latitude'],
                    'longitude': dict['gps_longitude'],
                    'altitude': dict['gps_altitude']}

    def config(self, date=None):
        """Get station config

        Retrieve either the latest, or a config for a specific date.

        :param date: date object for which to get the config
        :return: the full config for the station

        """
        if date is None:
            date = datetime.date.today()
        path = (self.urls['configuration']
                .format(station_number=self.station,
                        year=date.year, month=date.month, day=date.day))

        return self._get_json(path)

    def n_events(self, year='', month='', day='', hour=''):
        """Get number of events

        Note that it is possible to give only the year to get the total
        number of events in that year. If both year and month are given,
        the total events in that month are returned.

        :param year,month,day,hour: the date and time for which to
            get the number. It is possible to be less specific.
        :return: the number of events recorded by the station on date.

        """
        if year == '' and (month != '' or day != '' or hour != ''):
            raise Exception('You must also specify the year')
        elif month == '' and (day != '' or hour != ''):
            raise Exception('You must also specify the month')
        elif day == '' and hour != '':
            raise Exception('You must also specify the day')

        path = (self.urls['number_of_events']
                .format(station_number=self.station, year=year, month=month,
                        day=day, hour=hour).strip("/"))
        return self._get_json(path)

    def has_data(self, year='', month='', day=''):
        """Check for HiSPARC data

        :param year,month,day: the date for which to check. It is
            possible to be less specific.
        :return: boolean, indicating wether the station had air shower
            data on the date.

        """
        if year == '' and (month != '' or day != ''):
            raise Exception('You must also specify the year')
        elif month == '' and day != '':
            raise Exception('You must also specify the month')

        path = (self.urls['has_data'].format(station_number=self.station,
                                             year=year, month=month, day=day)
                .strip("/"))
        return self._get_json(path)

    def has_weather(self, year='', month='', day=''):
        """Check for weather data

        :param year,month,day: the date for which to check. It is
            possible to be less specific.
        :return: boolean, indicating wether the station had weather data
            on the date.

        """
        if year == '' and (month != '' or day != ''):
            raise Exception('You must also specify the year')
        elif month == '' and day != '':
            raise Exception('You must also specify the month')

        path = (self.urls['has_weather']
                .format(station_number=self.station,
                        year=year, month=month, day=day).strip("/"))
        return self._get_json(path)

    def event_trace(self, timestamp, nanoseconds):
        """Get the traces for a specific event

        The exact timestamp and nanoseconds for the event have to be
        given.

        :param timestamp,nanoseconds: the extended timestamp for which
            to get the traces
        :return: an array with the traces for each detector in ADCcounts

        """
        ext_timestamp = '%d%09d' % (timestamp, nanoseconds)
        path = (self.urls['event_trace'].format(station_number=self.station,
                                                ext_timestamp=ext_timestamp)
                .strip("/"))
        return self._get_json(path)

    def event_time(self, year, month, day):
        """Get the number of events per hour histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('hour', 'counts')
        path = self.src_urls['eventtime'].format(station_number=self.station,
                                                 year=year, month=month,
                                                 day=day)
        return self._get_csv(path, names=columns)

    def pulse_height(self, year, month, day):
        """Get the pulseheight histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('pulseheight', 'ph1', 'ph2', 'ph3', 'ph4')
        path = self.src_urls['pulseheight'].format(station_number=self.station,
                                                   year=year, month=month,
                                                   day=day)
        return self._get_csv(path, names=columns)

    def pulse_integral(self, year, month, day):
        """Get the pulseintegral histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('pulseintegral', 'pi1', 'pi2', 'pi3', 'pi4')
        path = self.src_urls['integral'].format(station_number=self.station,
                                                year=year, month=month,
                                                day=day)
        return self._get_csv(path, names=columns)

    def barometer(self, year, month, day):
        """Get the barometer dataset

        :param year,month,day: the date for which to get the dataset.
        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'air_pressure')
        path = self.src_urls['barometer'].format(station_number=self.station,
                                                 year=year, month=month,
                                                 day=day)
        return self._get_csv(path, names=columns)

    def temperature(self, year, month, day):
        """Get the temperature dataset

        :param year,month,day: the date for which to get the dataset.
        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'temperature')
        path = self.src_urls['temperature'].format(station_number=self.station,
                                                   year=year, month=month,
                                                   day=day)
        return self._get_csv(path, names=columns)

    @lazy
    def voltages(self):
        """Get the PMT voltage data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'voltage1', 'voltage2', 'voltage3', 'voltage4')
        path = self.src_urls['voltage'].format(station_number=self.station)
        return self._get_csv(path, names=columns)

    def voltage(self, timestamp):
        """Get PMT coltage data for specific timestamp

        :param timestamp: timestamp for which the value is valid.
        :return: list of values for given timestamp.

        """
        voltages = self.voltages
        idx = self.get_active_index(voltages['timestamp'], timestamp)
        voltage = (voltages[idx]['voltage1'], voltages[idx]['voltage2'],
                   voltages[idx]['voltage3'], voltages[idx]['voltage4'])
        return voltage

    @lazy
    def currents(self):
        """Get the PMT current data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'current1', 'current2', 'current3', 'current4')
        path = self.src_urls['current'].format(station_number=self.station)
        return self._get_csv(path, names=columns)

    def current(self, timestamp):
        """Get PMT current data for specific timestamp

        :param timestamp: timestamp for which the value is valid.
        :return: list of values for given timestamp.

        """
        currents = self.currents
        idx = self.get_active_index(currents['timestamp'], timestamp)
        current = (currents[idx]['current1'], currents[idx]['current2'],
                   currents[idx]['current3'], currents[idx]['current4'])
        return current

    @lazy
    def gps_locations(self):
        """Get the GPS location data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'latitude', 'longitude', 'altitude')
        path = self.src_urls['gps'].format(station_number=self.station)
        return self._get_csv(path, names=columns)

    def gps_location(self, timestamp):
        """Get GPS location for specific timestamp

        :param timestamp: timestamp for which the value is valid.
        :return: list of values for given timestamp.

        """
        locations = self.gps_locations
        idx = self.get_active_index(locations['timestamp'], timestamp)
        location = {'latitude': locations[idx]['latitude'],
                    'longitude': locations[idx]['longitude'],
                    'altitude': locations[idx]['altitude']}
        return location
