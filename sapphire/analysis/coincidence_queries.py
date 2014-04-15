import itertools
import warnings

import tables

from sapphire import api


class CoincidenceQuery(object):

    """Perform queries on an ESD file where coincidences have been analysed.

    Functions in this class build and perform queries to easily filter
    coincidences based on station participation and station
    organization. This assumes coincidences have been analysed using
    :class:`sapphire.analysis.CoincidenceESD`.

    First initiate the class so it can get the correct tables, then get
    coincidences using one of the functions (all_coincidences, any, all,
    at_least), then get the events for these coincidences using other
    functions (all_events, events_from_[stations, subcluster, cluster]).

    An exception will occur when you include a station in a query that
    does not occur in the datafile.

    Example usage::

        from sapphire.analysis.coincidence_queries import CoincidenceQuery

        cq = CoincidenceQuery('2013_8_1.h5')
        coincidences = cq.all([501, 502, 503, 504])
        events = cq.all_events(coincidences)
        specific_events = cq.events_from_stations(coincidences,
                                                  [501, 502, 503, 504])
        cq.data.close()

    """

    def __init__(self, data_path, coincidence_group='/coincidences'):
        self.data = tables.openFile(data_path, 'r')
        self.coincidences = self.data.getNode(coincidence_group,
                                              'coincidences')
        self.s_index = self.data.getNode(coincidence_group, 's_index')
        self.c_index = self.data.getNode(coincidence_group, 'c_index')

    def all_coincidences(self):
        """Get all coincidences

        :return: all coincidences.

        """
        coincidences = self.coincidences.read()
        return coincidences

    def any(self, stations):
        """Filter for coincidences that contain any of the given stations

        :param stations: list of stations from which any need to be in
                         a coincidence.
        :return: coincidences matching the query.

        """
        s_columns = self._get_s_columns(stations)
        any_query = ' | '.join(s_columns)
        filtered_coincidences = self.coincidences.readWhere(any_query)
        return filtered_coincidences

    def all(self, stations):
        """Filter for coincidences that contain all of the given stations

        :param stations: list of stations which all need to be in a
                         coincidence.
        :return: coincidences matching the query.

        """
        s_columns = self._get_s_columns(stations)
        all_query = ' & '.join(s_columns)
        filtered_coincidences = self.coincidences.readWhere(all_query)
        return filtered_coincidences

    def at_least(self, stations, n):
        """Filter coincidences to contain at least n of the given stations

        :param stations: list of stations from which any at least n of
                         any combinations need to be in a coincidence.
        :param n: minimum number of given stations to be in a coincidence.
        :return: coincidences matching the query.

        """
        s_columns = self._get_s_columns(stations)
        s_combinations = ['(%s)' % (' & '.join(combo))
                          for combo in itertools.combinations(s_columns, n)]
        at_least_query = ' | '.join(s_combinations)
        filtered_coincidences = self.coincidences.readWhere(at_least_query)
        return filtered_coincidences

    def _get_s_columns(self, stations):
        """Get column names for given stations

        :param stations: list of station numbers.
        :return: list of strings with column titles for each station.

        """
        return ['s%d' % station for station in stations]

    def _get_events(self, coincidence):
        """Get events belonging to a coincidence

        :param coincidence: A coincidence row.
        :return: list of tuples containing station numbers and events.

        """
        events = []
        c_idx = self.c_index[coincidence['id']]
        for s_idx, e_idx in c_idx:
            s_path = self.s_index[s_idx]
            station_number = int(s_path.split('station_')[-1])
            s_group = self.data.getNode(s_path)
            events.append((station_number, s_group.events[e_idx]))
        return events

    def all_events(self, coincidences):
        """Get events for all coincidences.

        :param coincidences: list of coincidence rows.
        :return: list of events for each coincidence.

        """
        coincidence_events = [self._get_events(coincidence)
                              for coincidence in coincidences]
        return coincidence_events

    def minimum_events_for_coincidence(self, coincidences_events, n=2):
        """Filter coincidences to only include those with at least n events.

        :param coincidences_events: list of events for each coincidence.
        :param n: minimum number of events per coincidence.

        """
        filtered_coincidences = [coincidence
                                 for coincidence in coincidences_events
                                 if len(coincidence) >= n]
        return filtered_coincidences

    def events_from_stations(self, coincidences, stations):
        """Only get events for specific stations for coincidences.

        :param coincidences: list of coincidence rows.
        :param stations: list of station numbers to filter events for.
        :return: list of filtered events for each coincidence.

        """
        events_iterator = (self._get_events(coincidence)
                           for coincidence in coincidences)
        coincidences_events = [self._events_from_stations(events, stations)
                               for events in events_iterator]
        return self.minimum_events_for_coincidence(coincidences_events)

    def _events_from_stations(self, events, stations):
        """Get only events from the chosen stations

        :param events: list of tuples containing a station number and an
                       event.
        :param stations: list of station numbers to filter events for.
        :return: tuples of station numbers and events for matched events.

        """
        filtered_events = [event for event in events if event[0] in stations]
        return filtered_events

    def events_in_subcluster(self, coincidences, subcluster):
        """Filter to only events from stations in a specific subcluster

        :param coincidences: list of events for each coincidence.
        :param subcluster: subcluster number to filter events for.
        :return: tuples of station numbers and events for matched events.

        """
        network = api.Network()
        stations = network.station_numbers(subcluster=subcluster)
        filtered_events = self.events_from_stations(coincidences)

        return filtered_events

    def events_in_cluster(self, coincidences, cluster):
        """Filter to only events from stations in a specific subcluster

        :param coincidences: list of events for each coincidence.
        :param cluster: cluster number to filter events for.
        :return: tuples of station numbers and events for matched events.

        """
        network = api.Network()
        stations = network.stations_numbers(cluster=cluster)
        filtered_events = self.events_from_stations(coincidences, stations)

        return filtered_events
