import re
import itertools
import warnings

import tables

from .. import api


class CoincidenceQuery(object):

    """Perform queries on an ESD file where coincidences have been analysed.

    Functions in this class build and perform queries to easily filter
    coincidences based on station participation and station
    organization. This assumes coincidences have been analysed using
    :class:`~sapphire.analysis.coincidences.CoincidencesESD`.

    First initiate the class so it can get the correct tables, then get
    coincidences using one of the functions (all_coincidences, any, all,
    at_least), then get the events for these coincidences using other
    functions (all_events, events_from_[stations, subcluster, cluster]).

    An exception will occur when you include a station in a query that
    does not occur in the datafile.

    Example usage::

        >>> from sapphire import CoincidenceQuery
        >>> cq = CoincidenceQuery('2013_8_1.h5')
        >>> coincidences = cq.all([501, 502, 503, 504], iterator=True)
        >>> events = cq.all_events(coincidences)
        >>> specific_events = cq.events_from_stations(coincidences,
        ...                                          [501, 502, 503, 504])
        >>> cq.finish()

    """

    def __init__(self, data, coincidence_group='/coincidences'):
        """Setup variables to point to the tables

        :param data: either a PyTables file or path to a HDF5 file.
        :param coincidence_group: path to the coincidences group.

        """
        if not isinstance(data, tables.File):
            self.data = tables.open_file(data, 'r')
        else:
            self.data = data
        self.coincidences = self.data.get_node(coincidence_group,
                                               'coincidences')
        self.c_index = self.data.get_node(coincidence_group, 'c_index')
        self.s_index = self.data.get_node(coincidence_group, 's_index')
        self.s_nodes = []
        for s_path in self.s_index:
            try:
                self.s_nodes.append(self.data.get_node(s_path))
            except tables.NoSuchNodeError:
                warnings.warn('Missing some station groups. This is no '
                              'problem if those are not in coincidences.')
                self.s_nodes.append(None)
        re_number = re.compile('[0-9]+$')
        self.s_numbers = [int(re_number.search(s_path).group())
                          for s_path in self.s_index]

        try:
            self.reconstructions = self.data.get_node(coincidence_group,
                                                      'reconstructions')
            self.reconstructed = True
        except tables.NoSuchNodeError:
            self.reconstructed = False

    def finish(self):
        """Clean-up after using

        Do not use if you opened the file yourself and still intent to
        use it.

        """
        self.data.close()

    def all_coincidences(self, iterator=False):
        """Get all coincidences

        :return: all coincidences.

        """
        if iterator:
            coincidences = self.coincidences.iterrows()
        else:
            coincidences = self.coincidences.read()
        return coincidences

    def any(self, stations, start=None, stop=None, iterator=False):
        """Filter for coincidences that contain any of the given stations

        :param stations: list of stations from which any need to be in
                         a coincidence.
        :return: coincidences matching the query.

        """
        s_columns = self._get_allowed_s_columns(stations)
        if len(s_columns) == 0:
            # no stations would result in a bad query string
            return []
        query = '(%s)' % ' | '.join(s_columns)
        query = self._add_timestamp_filter(query, start, stop)
        filtered_coincidences = self.perform_query(query, iterator)
        return filtered_coincidences

    def all(self, stations, start=None, stop=None, iterator=False):
        """Filter for coincidences that contain all of the given stations

        :param stations: list of stations which all need to be in a
                         coincidence.
        :return: coincidences matching the query.

        """
        s_columns = self._get_allowed_s_columns(stations)
        if len(s_columns) < len(stations):
            # Not all requested stations exist in the data so it's
            # impossible to find a coincidence with all stations.
            return []
        query = '(%s)' % ' & '.join(s_columns)
        query = self._add_timestamp_filter(query, start, stop)
        filtered_coincidences = self.perform_query(query, iterator)
        return filtered_coincidences

    def at_least(self, stations, n, start=None, stop=None, iterator=False):
        """Filter coincidences to contain at least n of the given stations

        :param stations: list of stations from which any at least n of
                         any combinations need to be in a coincidence.
        :param n: minimum number of given stations to be in a coincidence.
        :return: coincidences matching the query.

        """
        s_columns = self._get_allowed_s_columns(stations)
        if len(s_columns) < n:
            # No combinations possible because there are to few stations
            return []
        s_combinations = ['(%s)' % (' & '.join(combo))
                          for combo in itertools.combinations(s_columns, n)]
        query = '(%s)' % ' | '.join(s_combinations)
        query = self._add_timestamp_filter(query, start, stop)
        filtered_coincidences = self.perform_query(query, iterator)
        return filtered_coincidences

    def timerange(self, start, stop, iterator=False):
        """Query based on timestamps

        :param start: timestamp from which to look for coincidences.
        :param stop: end timestamp for coincidences.
        :return: coincidences within the specified timerange.

        """
        query = '(%d <= timestamp) & (timestamp < %d)' % (start, stop)
        filtered_coincidences = self.perform_query(query, iterator)
        return filtered_coincidences

    def _add_timestamp_filter(self, query, start=None, stop=None):
        """Add timestamp filter to the query

        :param start: timestamp from which to look for coincidences.
        :param stop: end timestamp for coincidences.
        :return: query with added timestamp filter.

        """
        if start:
            query += ' & (%d <= timestamp)' % start
        if stop:
            query += ' & (timestamp < %d)' % stop

        return query

    def perform_query(self, query, iterator=False):
        """Perform a query on the coincidences table

        :param query: a valid PyTables query string for the coincidences
                      table.
        :return: coincidences matching the query.

        """
        if iterator:
            filtered_coincidences = self.coincidences.where(query)
        else:
            filtered_coincidences = self.coincidences.read_where(query)
        return filtered_coincidences

    def _get_allowed_s_columns(self, stations):
        """Get column names for given stations

        This ensures that the columnnames actually exist and are unique,
        otherwise an exception may be raised.

        :param stations: list of station numbers.
        :return: list of strings with column titles for each station,
                 but only if the column exists.

        """
        s_columns = self._get_s_columns(stations)
        return set(self.coincidences.colnames).intersection(s_columns)

    @staticmethod
    def _get_s_columns(stations):
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
            station_number = self.s_numbers[s_idx]
            s_node = self.s_nodes[s_idx]
            if s_node is None:
                warnings.warn('Missing station group for station id %d. '
                              'Events from it are excluded.' % s_idx)
                continue
            events.append((station_number, s_node.events[e_idx]))
        return events

    def _get_reconstructions(self, coincidence):
        """Get event reconstructions belonging to a coincidence

        :param coincidence: A coincidence row.
        :return: list of tuples containing station numbers and
                 reconstructed events.

        """
        reconstructions = []
        c_idx = self.c_index[coincidence['id']]
        for s_idx, e_idx in c_idx:
            station_number = self.s_numbers[s_idx]
            s_node = self.s_nodes[s_idx]
            if s_node is None:
                warnings.warn('Missing station group for station id %d.'
                              'Reconstructions from it are excluded.' % s_idx)
                continue
            rec_table = s_node.reconstructions
            reconstructions.append((station_number, rec_table[e_idx]))
        return reconstructions

    def _get_reconstruction(self, coincidence):
        """Get coincidence reconstruction belonging to a coincidence

        :param coincidence: A coincidence row.
        :return: reconstructed coincidence.

        """
        if self.reconstructed:
            reconstruction = self.reconstructions[coincidence['id']]
            return reconstruction
        else:
            raise Exception('Coincidences are not (properly) reconstructed.'
                            'Perform reconstructions and reinitialize this '
                            'class.')

    def all_events(self, coincidences, n=0):
        """Get all events for the given coincidences.

        :param coincidences: list of coincidence rows.
        :param n: minimum number of events per coincidence.
        :return: list of events for each coincidence.

        """
        coincidence_events = (self._get_events(coincidence)
                              for coincidence in coincidences)
        return self.minimum_events_for_coincidence(coincidence_events, n)

    def all_reconstructions(self, coincidences, n=0):
        """Get all reconstructed events for the given coincidences.

        :param coincidences: list of coincidence rows.
        :param n: minimum number of events per coincidence.
        :return: list of reconstructed events for each coincidence.

        """
        coincidence_recs = (self._get_reconstructions(coincidence)
                            for coincidence in coincidences)
        return self.minimum_events_for_coincidence(coincidence_recs, n)

    def minimum_events_for_coincidence(self, coincidences_events, n=2):
        """Filter coincidences to only include those with at least n events.

        :param coincidences_events: list of events for each coincidence.
        :param n: minimum number of events per coincidence.

        """
        filtered_coincidences = (coincidence
                                 for coincidence in coincidences_events
                                 if len(coincidence) >= n)
        return filtered_coincidences

    def events_from_stations(self, coincidences, stations, n=2):
        """Only get events for specific stations for coincidences.

        :param coincidences: list of coincidence rows.
        :param stations: list of station numbers to filter events for.
        :return: list of filtered events for each coincidence.

        """
        events_iterator = (self._get_events(coincidence)
                           for coincidence in coincidences)
        coincidences_events = (self._events_from_stations(events, stations)
                               for events in events_iterator)
        return self.minimum_events_for_coincidence(coincidences_events, n)

    def reconstructions_from_stations(self, coincidences, stations, n=2):
        """Only get reconstructions for specific stations for coincidences.

        :param coincidences: list of coincidence rows.
        :param stations: list of station numbers to filter events for.
        :return: list of filtered reconstructed events for each coincidence.

        """
        reconstructions_iterator = (self._get_reconstructions(coincidence)
                                    for coincidence in coincidences)
        coincidences_recs = (self._events_from_stations(recs, stations)
                             for recs in reconstructions_iterator)
        return self.minimum_events_for_coincidence(coincidences_recs, n)

    def _events_from_stations(self, events, stations):
        """Get only events from the chosen stations

        :param events: list of tuples containing a station number and an
                       event.
        :param stations: list of station numbers to filter events for.
        :return: tuples of station numbers and events for matched events.

        """
        filtered_events = [event for event in events if event[0] in stations]
        return filtered_events

    def events_in_subcluster(self, coincidences, subcluster, n=2):
        """Filter to only events from stations in a specific subcluster

        :param coincidences: list of events for each coincidence.
        :param subcluster: subcluster number to filter events for.
        :return: tuples of station numbers and events for matched events.

        """
        network = api.Network()
        stations = network.station_numbers(subcluster=subcluster)
        filtered_events = self.events_from_stations(coincidences, stations, n)

        return filtered_events

    def events_in_cluster(self, coincidences, cluster, n=2):
        """Filter to only events from stations in a specific subcluster

        :param coincidences: list of events for each coincidence.
        :param cluster: cluster number to filter events for.
        :return: tuples of station numbers and events for matched events.

        """
        network = api.Network()
        stations = network.station_numbers(cluster=cluster)
        filtered_events = self.events_from_stations(coincidences, stations, n)

        return filtered_events
