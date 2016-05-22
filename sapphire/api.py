""" Access the HiSPARC public database API.

    This provides easy classes and functions to access the HiSPARC
    publicdb API. This takes care of the url retrieval and conversion
    from JSON to Python dictionaries.

    Example usage:

    .. code-block:: python

        >>> from sapphire import Station
        >>> stations = [5, 3102, 504, 7101, 8008, 13005]
        >>> clusters = [Station(station).cluster for station in stations]
        >>> for station, cluster in zip(stations, clusters):
        ...     print 'Station %d is in cluster %s.' % (station, cluster)
        Station 5 is in cluster Amsterdam.
        Station 3102 is in cluster Leiden.
        Station 504 is in cluster Amsterdam.
        Station 7101 is in cluster Enschede.
        Station 8008 is in cluster Eindhoven.
        Station 13005 is in cluster Bristol.

"""
import logging
import datetime
import json
import warnings
from os import path, extsep
from urllib2 import urlopen, HTTPError, URLError
from StringIO import StringIO

from lazy import lazy
from numpy import (genfromtxt, atleast_1d, zeros, ones, logical_and,
                   count_nonzero, negative)

from .utils import get_active_index, memoize
from .transformations.clock import process_time

logger = logging.getLogger('api')

API_BASE = 'http://data.hisparc.nl/api/'
SRC_BASE = 'http://data.hisparc.nl/show/source/'
LOCAL_BASE = path.join(path.dirname(__file__), 'data')


class API(object):

    """Base API class

    This provided the methods to retrieve data from the API. The results
    are converted from JSON data to Python objects (dict/list/etc).
    Support is also provided for the retrieval of Source TSV data, which
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
        'azimuth': 'azimuth/{station_number}/{year}/{month}/{day}/',
        'zenith': 'zenith/{station_number}/{year}/{month}/{day}/',
        'barometer': 'barometer/{station_number}/{year}/{month}/{day}/',
        'temperature': 'temperature/{station_number}/{year}/{month}/{day}/',
        'electronics': 'electronics/{station_number}/',
        'voltage': 'voltage/{station_number}/',
        'current': 'current/{station_number}/',
        'gps': 'gps/{station_number}/',
        'trigger': 'trigger/{station_number}/',
        'layout': 'layout/{station_number}/',
        'detector_timing_offsets': 'detector_timing_offsets/{station_number}/',
        'station_timing_offsets': 'station_timing_offsets/{station_1}/'
                                  '{station_2}/'}

    def __init__(self, force_fresh=False, force_stale=False):
        """Initialize API class

        :param force_fresh,force_stale: if either of these is set to True the
            data must either loaded from server or from local data. Be default
            fresh data is prefered, but falls back to local data.

        """
        self.force_fresh = force_fresh
        self.force_stale = force_stale

    def _get_json(self, urlpath):
        """Retrieve a JSON from the HiSPARC API

        :param urlpath: api urlpath to retrieve (i.e. after API_BASE).
        :return: the data returned by the api as dictionary or integer.

        """
        if self.force_fresh and self.force_stale:
            raise Exception('Can not force fresh and stale simultaneously.')
        try:
            if self.force_stale:
                raise Exception
            json_data = self._retrieve_url(urlpath)
            data = json.loads(json_data)
        except Exception:
            if self.force_fresh:
                raise Exception('Couldn\'t get requested data from server.')
            localpath = path.join(LOCAL_BASE,
                                  urlpath.strip('/') + extsep + 'json')
            try:
                with open(localpath) as localdata:
                    data = json.load(localdata)
            except:
                if self.force_stale:
                    raise Exception('Couldn\'t find requested data locally.')
                raise Exception('Couldn\'t get requested data from server '
                                'nor find it locally.')
            if not self.force_stale:
                warnings.warn('Using local data. Possibly outdated.')

        return data

    def _get_tsv(self, urlpath, names=None):
        """Retrieve a Source TSV from the HiSPARC Public Database

        :param urlpath: tsv urlpath to retrieve (i.e. path after SRC_BASE).
        :param names: data column names.
        :return: the data returned as array.

        """
        if self.force_fresh and self.force_stale:
            raise Exception('Can not force fresh and stale simultaneously.')
        try:
            if self.force_stale:
                raise Exception
            tsv_data = self._retrieve_url(urlpath, base=SRC_BASE)
        except Exception:
            if self.force_fresh:
                raise Exception('Couldn\'t get requested data from server.')
            localpath = path.join(LOCAL_BASE,
                                  urlpath.strip('/') + extsep + 'tsv')
            try:
                with warnings.catch_warnings():
                    warnings.filterwarnings('ignore')
                    data = genfromtxt(localpath, delimiter='\t', dtype=None,
                                      names=names)
            except:
                if self.force_stale:
                    raise Exception('Couldn\'t find requested data locally.')
                raise Exception('Couldn\'t get requested data from server '
                                'nor find it locally.')
            if not self.force_stale:
                warnings.warn('Using local data. Possibly outdated.')
        else:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                data = genfromtxt(StringIO(tsv_data), delimiter='\t',
                                  dtype=None, names=names)

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
    def validate_partial_date(year='', month='', day='', hour=''):
        if year == '' and (month != '' or day != '' or hour != ''):
            raise Exception('You must also specify the year')
        elif month == '' and (day != '' or hour != ''):
            raise Exception('You must also specify the month')
        elif day == '' and hour != '':
            raise Exception('You must also specify the day')


class Network(API):

    """Get info about the network (countries/clusters/subclusters/stations)"""

    @lazy
    def _all_countries(self):
        """All countries data"""

        path = self.urls['countries']
        return self._get_json(path)

    def countries(self):
        """Get a list of countries

        :return: all countries in the region

        """
        return self._all_countries

    def country_numbers(self):
        """Same as countries but only retuns a list of country numbers"""

        countries = self.countries()
        return [country['number'] for country in countries]

    @lazy
    def _all_clusters(self):
        """All countries data"""

        path = self.urls['clusters']
        return self._get_json(path)

    def clusters(self, country=None):
        """Get a list of clusters

        :param country: the number of the country for which to get all
            clusters.
        :return: all clusters in the region

        """
        self.validate_numbers(country)
        if country is None:
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

    @lazy
    def _all_subclusters(self):
        """All countries data"""

        path = self.urls['subclusters']
        return self._get_json(path)

    def subclusters(self, country=None, cluster=None):
        """Get a list of subclusters

        :param country,cluster: the number of the region for which to get
            the subclusters it contains, only one or none should
            be specified.
        :return: all subclusters in the region

        """
        self.validate_numbers(country, cluster)
        if country is None and cluster is None:
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

    @lazy
    def _all_stations(self):
        """All stations data"""

        path = self.urls['stations']
        return self._get_json(path)

    def stations(self, country=None, cluster=None, subcluster=None):
        """Get a list of stations

        :param country,cluster,subcluster: the number of the region
            for which to get all stations, only one or none should
            be specified.
        :return: all stations in the region

        """
        self.validate_numbers(country, cluster, subcluster)
        if country is None and cluster is None and subcluster is None:
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
        elif subcluster is not None:
            path = (self.urls['stations_in_subcluster']
                    .format(subcluster_number=subcluster))
            stations = self._get_json(path)
        return stations

    def station_numbers(self, country=None, cluster=None, subcluster=None):
        """Same as stations but only retuns a list of station numbers"""

        stations = self.stations(country=country, cluster=cluster,
                                 subcluster=subcluster)
        return [station['number'] for station in stations]

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

    def stations_with_data(self, year='', month='', day=''):
        """Get a list of stations with data on the specified date

        :param year,month,day: the date for which to check. It is
            possible to be less specific.
        :return: all stations with data.

        """
        self.validate_partial_date(year, month, day)

        path = (self.urls['stations_with_data']
                .format(year=year, month=month, day=day).strip("/"))
        return self._get_json(path)

    def stations_with_weather(self, year='', month='', day=''):
        """Get a list of stations with weather data on the specified date

        :param year,month,day: the date for which to check. It is
            possible to be less specific.
        :return: all stations with weather data.

        """
        self.validate_partial_date(year, month, day)

        path = (self.urls['stations_with_weather']
                .format(year=year, month=month, day=day).strip("/"))
        return self._get_json(path)

    def coincidence_time(self, year, month, day):
        """Get the coincidences per hour histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('hour', 'counts')
        path = self.src_urls['coincidencetime'].format(year=year, month=month,
                                                       day=day)
        return self._get_tsv(path, names=columns)

    def coincidence_number(self, year, month, day):
        """Get the number of stations in coincidence histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('n', 'counts')
        path = self.src_urls['coincidencenumber'].format(year=year,
                                                         month=month,
                                                         day=day)
        return self._get_tsv(path, names=columns)

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

    def uptime(self, stations, start=None, end=None):
        """Get number of hours which stations have been simultaneously active

        Using hourly eventrate data the number of hours in which the given
        stations all had data is determined. Only hours in which each station
        had a reasonable eventrate are counted, to exclude bad data.

        :param stations: station number or a list of station numbers.
        :param start,end: start, end timestamp.
        :returns: number of hours with simultaneous data.

        """
        data = {}

        if not hasattr(stations, '__len__'):
            stations = [stations]

        for sn in stations:
            data[sn] = Station(sn, force_fresh=self.force_fresh,
                               force_stale=self.force_stale).event_time()

        first = min(values['timestamp'][0] for values in data.values())
        last = max(values['timestamp'][-1] for values in data.values())

        len_array = (last - first) / 3600 + 1
        all_active = ones(len_array)

        for sn in data.keys():
            is_active = zeros(len_array)
            start_i = (data[sn]['timestamp'][0] - first) / 3600
            end_i = start_i + len(data[sn])
            is_active[start_i:end_i] = (data[sn]['counts'] > 500) &\
                                       (data[sn]['counts'] < 5000)
            all_active = logical_and(all_active, is_active)

        # filter start, end
        if start is not None:
            start_index = max(0, process_time(start) - first) / 3600
        else:
            start_index = 0

        if end is not None:
            end_index = min(last, process_time(end) - first) / 3600
        else:
            end_index = len(all_active)

        return count_nonzero(all_active[start_index:end_index])


class Station(API):

    """Access data about a single station"""

    def __init__(self, station, date=None, force_fresh=False,
                 force_stale=False):
        """Initialize station

        :param station: station number.
        :param date: date object for which to get the station information.
        :param force_fresh: set to True to require data to be fresh
                            from the server.
        :param force_stale: set to True to require data to be taken from local
                            data, not valid for all methods.

        """
        if force_fresh and force_stale:
            raise Exception('Can not force fresh and stale simultaneously.')
        if station not in Network(force_fresh=force_fresh,
                                  force_stale=force_stale).station_numbers():
            warnings.warn('Possibly invalid station, or without config.')
        self.force_fresh = force_fresh
        self.force_stale = force_stale
        self.station = station
        if date is None:
            self.year, self.month, self.day = ('', '', '')
        else:
            self.year, self.month, self.day = (date.year, date.month, date.day)

    @lazy
    def info(self):
        """Get general station info"""

        path = (self.urls['station_info']
                .format(station_number=self.station, year=self.year,
                        month=self.month, day=self.day).strip("/"))
        return self._get_json(path)

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
            station = Station(self.station, date, self.force_fresh,
                              self.force_stale)
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
        self.validate_partial_date(year, month, day, hour)

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
        self.validate_partial_date(year, month, day)

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
        self.validate_partial_date(year, month, day)

        path = (self.urls['has_weather']
                .format(station_number=self.station,
                        year=year, month=month, day=day).strip("/"))
        return self._get_json(path)

    def event_trace(self, timestamp, nanoseconds, raw=False):
        """Get the traces for a specific event

        The exact timestamp and nanoseconds for the event have to be
        given.

        :param timestamp,nanoseconds: the extended timestamp for which
            to get the traces.
        :param raw: get the raw trace, without the subtracted baselines.
        :return: an array with the traces for each detector in ADCcounts

        """
        ext_timestamp = '%d%09d' % (timestamp, nanoseconds)
        path = (self.urls['event_trace'].format(station_number=self.station,
                                                ext_timestamp=ext_timestamp)
                .strip("/"))
        if raw is True:
            path += '?raw'
        return self._get_json(path)

    def event_time(self, year='', month='', day=''):
        """Get the number of events per hour histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('timestamp', 'counts')
        path = self.src_urls['eventtime'].format(station_number=self.station,
                                                 year=year, month=month,
                                                 day=day).strip("/")
        return self._get_tsv(path, names=columns)

    def pulse_height(self, year, month, day):
        """Get the pulseheight histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('pulseheight', 'ph1', 'ph2', 'ph3', 'ph4')
        path = self.src_urls['pulseheight'].format(station_number=self.station,
                                                   year=year, month=month,
                                                   day=day)
        return self._get_tsv(path, names=columns)

    def pulse_integral(self, year, month, day):
        """Get the pulseintegral histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('pulseintegral', 'pi1', 'pi2', 'pi3', 'pi4')
        path = self.src_urls['integral'].format(station_number=self.station,
                                                year=year, month=month,
                                                day=day)
        return self._get_tsv(path, names=columns)

    def azimuth(self, year, month, day):
        """Get the azimuth histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('angle', 'counts')
        path = self.src_urls['azimuth'].format(station_number=self.station,
                                               year=year, month=month, day=day)
        return self._get_tsv(path, names=columns)

    def zenith(self, year, month, day):
        """Get the zenith histogram

        :param year,month,day: the date for which to get the histogram.
        :return: array of bins and counts.

        """
        columns = ('angle', 'counts')
        path = self.src_urls['zenith'].format(station_number=self.station,
                                              year=year, month=month, day=day)
        return self._get_tsv(path, names=columns)

    def barometer(self, year, month, day):
        """Get the barometer dataset

        :param year,month,day: the date for which to get the dataset.
        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'air_pressure')
        path = self.src_urls['barometer'].format(station_number=self.station,
                                                 year=year, month=month,
                                                 day=day)
        return self._get_tsv(path, names=columns)

    def temperature(self, year, month, day):
        """Get the temperature dataset

        :param year,month,day: the date for which to get the dataset.
        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'temperature')
        path = self.src_urls['temperature'].format(station_number=self.station,
                                                   year=year, month=month,
                                                   day=day)
        return self._get_tsv(path, names=columns)

    @lazy
    def electronics(self):
        """Get the electronics version data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'master', 'slave', 'master_fpga', 'slave_fpga')
        path = self.src_urls['electronics'].format(station_number=self.station)
        return self._get_tsv(path, names=columns)

    def electronic(self, timestamp):
        """Get electronics version data for specific timestamp

        :param timestamp: timestamp for which the values are valid.
        :return: list of values for given timestamp.

        """
        electronics = self.electronics
        idx = get_active_index(electronics['timestamp'], timestamp)
        electronic = [electronics[idx][field] for field in
                      ('master', 'slave', 'master_fpga', 'slave_fpga')]
        return electronic

    @lazy
    def voltages(self):
        """Get the PMT voltage data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'voltage1', 'voltage2', 'voltage3', 'voltage4')
        path = self.src_urls['voltage'].format(station_number=self.station)
        return self._get_tsv(path, names=columns)

    def voltage(self, timestamp):
        """Get PMT voltage data for specific timestamp

        :param timestamp: timestamp for which the values are valid.
        :return: list of values for given timestamp.

        """
        voltages = self.voltages
        idx = get_active_index(voltages['timestamp'], timestamp)
        voltage = [voltages[idx]['voltage%d' % i] for i in range(1, 5)]
        return voltage

    @lazy
    def currents(self):
        """Get the PMT current data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'current1', 'current2', 'current3', 'current4')
        path = self.src_urls['current'].format(station_number=self.station)
        return self._get_tsv(path, names=columns)

    def current(self, timestamp):
        """Get PMT current data for specific timestamp

        :param timestamp: timestamp for which the values are valid.
        :return: list of values for given timestamp.

        """
        currents = self.currents
        idx = get_active_index(currents['timestamp'], timestamp)
        current = [currents[idx]['current%d' % i] for i in range(1, 5)]
        return current

    @lazy
    def gps_locations(self):
        """Get the GPS location data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'latitude', 'longitude', 'altitude')
        path = self.src_urls['gps'].format(station_number=self.station)
        return self._get_tsv(path, names=columns)

    def gps_location(self, timestamp=None):
        """Get GPS location for specific timestamp

        :param timestamp: optional timestamp or datetime object for which
            the values are valid.
        :return: dictionary with the values for given timestamp.

        """
        if timestamp is None:
            timestamp = process_time(self.date)
        else:
            timestamp = process_time(timestamp)
        locations = self.gps_locations
        idx = get_active_index(locations['timestamp'], timestamp)
        location = {'latitude': locations[idx]['latitude'],
                    'longitude': locations[idx]['longitude'],
                    'altitude': locations[idx]['altitude']}
        return location

    @lazy
    def triggers(self):
        """Get the trigger config data

        :return: array of timestamps and values.

        """
        columns = ('timestamp',
                   'low1', 'low2', 'low3', 'low4',
                   'high1', 'high2', 'high3', 'high4',
                   'n_low', 'n_high', 'and_or', 'external')
        path = self.src_urls['trigger'].format(station_number=self.station)
        return self._get_tsv(path, names=columns)

    def trigger(self, timestamp):
        """Get trigger config for specific timestamp

        :param timestamp: timestamp for which the values are valid.
        :return: thresholds and trigger values for given timestamp.

        """
        triggers = self.triggers
        idx = get_active_index(triggers['timestamp'], timestamp)
        thresholds = [[triggers[idx]['%s%d' % (t, i)]
                       for t in ('low', 'high')]
                      for i in range(1, 5)]
        trigger = [triggers[idx][t]
                   for t in 'n_low', 'n_high', 'and_or', 'external']
        return thresholds, trigger

    @lazy
    def station_layouts(self):
        """Get the station layout data

        :return: array of timestamps and values.

        """
        columns = ('timestamp',
                   'radius1', 'alpha1', 'height1', 'beta1',
                   'radius2', 'alpha2', 'height2', 'beta2',
                   'radius3', 'alpha3', 'height3', 'beta3',
                   'radius4', 'alpha4', 'height4', 'beta4')
        base = self.src_urls['layout']
        path = base.format(station_number=self.station)
        return self._get_tsv(path, names=columns)

    def station_layout(self, timestamp):
        """Get station layout data for specific timestamp

        :param timestamp: timestamp for which the values are valid.
        :return: list of coordinates for given timestamp.

        """
        station_layouts = self.station_layouts
        idx = get_active_index(station_layouts['timestamp'], timestamp)
        station_layout = [[station_layouts[idx]['%s%d' % (c, i)]
                           for c in ('radius', 'alpha', 'height', 'beta')]
                          for i in range(1, 5)]
        return station_layout

    @lazy
    def detector_timing_offsets(self):
        """Get the detector timing offsets data

        :return: array of timestamps and values.

        """
        columns = ('timestamp', 'offset1', 'offset2', 'offset3', 'offset4')
        base = self.src_urls['detector_timing_offsets']
        path = base.format(station_number=self.station)
        return self._get_tsv(path, names=columns)

    def detector_timing_offset(self, timestamp):
        """Get detector timing offset data for specific timestamp

        :param timestamp: timestamp for which the values are valid.
        :return: list of values for given timestamp.

        """
        detector_timing_offsets = self.detector_timing_offsets
        idx = get_active_index(detector_timing_offsets['timestamp'],
                               timestamp)
        detector_timing_offset = [detector_timing_offsets[idx]['offset%d' % i]
                                  for i in range(1, 5)]

        return detector_timing_offset

    @memoize
    def station_timing_offsets(self, reference_station):
        """Get the station timing offset relative to reference_station

        :param reference_station: reference station
        :return: array of timestamps and values.

        """
        if reference_station == self.station:
            raise Exception('Reference station cannot be the same station')
        if reference_station > self.station:
            station_1, station_2 = self.station, reference_station
            toggle_sign = True
        else:
            station_2, station_1 = self.station, reference_station
            toggle_sign = False

        columns = ('timestamp', 'offset', 'rchi2')
        base = self.src_urls['station_timing_offsets']
        path = base.format(station_1=station_1, station_2=station_2)
        data = self._get_tsv(path, names=columns)
        if toggle_sign:
            data['offset'] = negative(data['offset'])
        return data

    def station_timing_offset(self, timestamp, reference_station):
        """Get station timing offset data for specific timestamp

        :param timestamp: timestamp for which the value is valid.
        :param reference_station: reference station
        :return: list of values for given timestamp.

        """
        station_timing_offsets = self.station_timing_offsets(reference_station)
        idx = get_active_index(station_timing_offsets['timestamp'],
                               timestamp)
        station_timing_offset = (station_timing_offsets[idx]['offset'],
                                 station_timing_offsets[idx]['rchi2'])

        return station_timing_offset
