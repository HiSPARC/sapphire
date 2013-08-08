""" HiSPARC api interface

    This provides easy classes and functions to access the HiSPARC
    publicdb API. This takes care of the url retrieval and conversion
    from JSON to Python dictionaries.

"""
import datetime
from urllib2 import urlopen, HTTPError, URLError
import json


API = {"stations": 'stations/',
       "stations_in_subcluster": 'subclusters/{subcluster_id}/',
       "subclusters": 'subclusters/',
       "subclusters_in_cluster": 'clusters/{cluster_id}/',
       "clusters": 'clusters/',
       "clusters_in_country": 'countries/{country_id}/',
       "countries": 'countries/',
       "stations_with_data": 'stations/data/{year}/{month}/{day}/',
       "stations_with_weather": 'stations/weather/{year}/{month}/{day}/',
       "station_info": 'station/{station_id}/',
       "has_data": 'station/{station_id}/data/{year}/{month}/{day}/',
       "has_weather": 'station/{station_id}/weather/{year}/{month}/{day}/',
       "configuration": 'station/{station_id}/config/{year}/{month}/{day}/',
       "number_of_events": 'station/{station_id}/num_events/{year}/{month}/{day}/{hour}/'}

__api = 'http://data.hisparc.nl/api/'


class Network(object):
    """Get info about the network (countries/clusters/subclusters/stations)"""

    def __init__(self):
        """Initialize network
        
        """
        path = API['countries']
        self.all_countries = _get_json(path)
        path = API['clusters']
        self.all_clusters = _get_json(path)
        path = API['subclusters']
        self.all_subclusters = _get_json(path)
        path = API['stations']
        self.all_stations = _get_json(path)

    def clusters(self, country=None):
        if country is None:
            return self.all_clusters
        else:
            stations = []
            path = API['clusters_in_country'].format(country_id=country)
            return _get_json(path)

    def subcluster(self, country=None, cluster=None):
        if country is None and cluster is None:
            return self.all_subclusters
        elif not country is None:
            subclusters = []
            path = API['clusters_in_country'].format(country_id=country)
            clusters = _get_json(path)
            for cluster in clusters:
                path = (API['subclusters_in_cluster']
                        .format(cluster_id=cluster['number']))
                subclusters.extend(_get_json(path))
            return subclusters
        else:
            path = API['subclusters_in_cluster'].format(cluster_id=cluster)
            return _get_json(path)


    def stations(self, country=None, cluster=None, subcluster=None):
        if country is None and cluster is None and subcluster is None:
            return self.all_stations
        elif not country is None:
            stations = []
            path = API['clusters_in_country'].format(country_id=country)
            clusters = _get_json(path)
            for cluster in clusters:
                path = (API['subclusters_in_cluster']
                        .format(cluster_id=cluster['number']))
                subclusters = _get_json(path)
                for subcluster in subclusters:
                    path = (API['stations_in_subcluster']
                            .format(subcluster_id=subcluster['number']))
                    stations.extend(_get_json(path))
            return stations
        elif not cluster is None:
            stations = []
            path = API['subclusters_in_cluster'].format(cluster_id=cluster)
            subclusters = _get_json(path)
            for subcluster in subclusters:
                path = (API['stations_in_subcluster']
                        .format(subcluster_id=subcluster['number']))
                stations.extend(_get_json(path))
            return stations
        elif not subcluster is None:
            path = (API['stations_in_subcluster']
                    .format(subcluster_id=subcluster))
            return _get_json(path)

class Station(object):
    """Access data about a single station"""

    def __init__(self, station):
        """Initialize station

        :param station: station number

        """
        self.station = station
        path = API['station_info'].format(station_id=self.station)
        self.info = _get_json(path)

    @property
    def country(self):
        return self.info['country']

    @property
    def cluster(self):
        return self.info['cluster']

    @property
    def subcluster(self):
        return self.info['subcluster']

    @property
    def n_detectors(self):
        return self.info['scintillators']

    def detectors(self, date=None):
        if date is None:
            dict = self.info
        else:
            raise Exception('Not supported yet')

        if dict['scintillators'] == 4:
            return {'scintillator1': dict['scintillator1'],
                    'scintillator2': dict['scintillator2'],
                    'scintillator3': dict['scintillator3'],
                    'scintillator4': dict['scintillator4']}
        else:
            return {'scintillator1': dict['scintillator1'],
                    'scintillator2': dict['scintillator2']}

    def location(self, date=None):
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
        if date is None:
            date = datetime.date.today()
        path = API['configuration'].format(station_id=self.station,
                                           year=date.year,
                                           month=date.month,
                                           day=date.day)

        return _get_json(path)


    def n_events(self, year='', month='', day='', hour=''):
        if year == '' and (month != '' or day != '' or hour != ''):
            raise Exception('You must also specify the year')
        elif month == '' and (day != '' or hour != ''):
            raise Exception('You must also specify the month')
        elif day == '' and hour != '':
            raise Exception('You must also specify the day')

        path = API['number_of_events'].format(station_id=self.station,
                                              year=year, month=month, day=day,
                                              hour=hour).strip("/")
        return _get_json(path)

    def has_data(self, year='', month='', day=''):
        if year == '' and (month != '' or day != ''):
            raise Exception('You must also specify the year')
        elif month == '' and day != '':
            raise Exception('You must also specify the month')

        path = (API['has_data'].format(station_id=self.station,
                                               year=year, month=month, day=day)
                                       .strip("/"))
        return _get_json(path)

    def has_weather(self, year='', month='', day=''):
        if year == '' and (month != '' or day != ''):
            raise Exception('You must also specify the year')
        elif month == '' and day != '':
            raise Exception('You must also specify the month')

        path = (API['has_weather'].format(station_id=self.station,
                                               year=year, month=month, day=day)
                                       .strip("/"))
        return _get_json(path)


def _get_json(urlpath):
    json_data = _retrieve_url(urlpath)
    data = json.loads(json_data)
    
    return data


def _retrieve_url(urlpath):

    url = __api + urlpath

    try:
        result = urlopen(url).read()
    except HTTPError, e:
        raise Exception('A HTTP %d error occured.' % e.code)
    except URLError:
        raise Exception('An URL error occured.')

    return result