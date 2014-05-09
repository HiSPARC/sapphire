""" Process HiSPARC events

    This module can be used analyse data to get observables like arrival
    times and particle count in each detector for each event.

    Example usage::

        import datetime

        import tables

        from sapphire.publicdb import download_data
        from sapphire.analysis import process_events


        STATIONS = [501, 503, 506]
        START = datetime.datetime(2013, 1, 1)
        END = datetime.datetime(2013, 1, 2)


        if __name__ == '__main__':
            station_groups = ['/s%d' % u for u in STATIONS]

            data = tables.open_file('data.h5', 'w')
            for station, group in zip(STATIONS, station_groups):
                download_data(data, group, station, START, END, get_blobs=True)
                proc = process_events.ProcessEvents(data, group)
                proc.process_and_store_results()
            data.close()

"""
import zlib
from itertools import izip
import operator

import tables
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
import progressbar as pb

from sapphire.analysis.find_mpv import FindMostProbableValueInSpectrum


ADC_THRESHOLD = 20
ADC_LOW_THRESHOLD = 253
ADC_HIGH_THRESHOLD = 323
ADC_TIME_PER_SAMPLE = 2.5e-9


class ProcessEvents(object):

    """Process HiSPARC events to obtain several observables.

    This class can be used to process a set of HiSPARC events and adds a
    few observables like particle arrival time and number of particles in
    the detector to a copy of the event table.

    """

    processed_events_description = {
        'event_id': tables.UInt32Col(pos=0),
        'timestamp': tables.Time32Col(pos=1),
        'nanoseconds': tables.UInt32Col(pos=2),
        'ext_timestamp': tables.UInt64Col(pos=3),
        'data_reduction': tables.BoolCol(pos=4),
        'trigger_pattern': tables.UInt32Col(pos=5),
        'baseline': tables.Int16Col(pos=6, shape=4, dflt=-1),
        'std_dev': tables.Int16Col(pos=7, shape=4, dflt=-1),
        'n_peaks': tables.Int16Col(pos=8, shape=4, dflt=-1),
        'pulseheights': tables.Int16Col(pos=9, shape=4, dflt=-1),
        'integrals': tables.Int32Col(pos=10, shape=4, dflt=-1),
        'traces': tables.Int32Col(pos=11, shape=4, dflt=-1),
        'event_rate': tables.Float32Col(pos=12),
        't1': tables.Float32Col(pos=13, dflt=-1),
        't2': tables.Float32Col(pos=14, dflt=-1),
        't3': tables.Float32Col(pos=15, dflt=-1),
        't4': tables.Float32Col(pos=16, dflt=-1),
        'n1': tables.Float32Col(pos=17, dflt=-1),
        'n2': tables.Float32Col(pos=18, dflt=-1),
        'n3': tables.Float32Col(pos=19, dflt=-1),
        'n4': tables.Float32Col(pos=20, dflt=-1),
        't_trigger': tables.Float32Col(pos=21, dflt=-1)}

    def __init__(self, data, group, source=None):
        """Initialize the class.

        :param data: the PyTables datafile
        :param group: the group containing the station data.  In normal
            cases, this is simply the group containing the events table.
        :param source: the name of the events table.  Default: None,
            meaning the default name 'events'.

        """
        self.data = data
        self.group = data.get_node(group)
        self.source = self._get_source(source)

    def process_and_store_results(self, destination=None, overwrite=False,
                                  limit=None):
        """Process events and store the results.

        :param destination: name of the table where the results will be
            written.  The default, None, corresponds to 'events'.
        :param overwrite: if True, overwrite previously obtained results.
        :param limit: the maximum number of events that will be stored.
            The default, None, corresponds to no limit.

        """
        self.limit = limit

        self._check_destination(destination, overwrite)

        self._clean_events_table()

        self._create_results_table()
        self._store_results_from_traces()
        self._store_number_of_particles()
        self._move_results_table_into_destination()

    def get_traces_for_event(self, event):
        """Return the traces from an event.

        :param event: a row from the events table.
        :returns: the traces: an array of pulseheight values.

        """
        traces = [self._get_trace(idx) for idx in event['traces'] if idx >= 0]

        # Make traces follow NumPy conventions
        traces = np.array(traces).T
        return traces

    def get_traces_for_event_index(self, idx):
        """Return the traces from event #idx.

        :param idx: the index number of the event.
        :returns: the traces: an array of pulseheight values.

        """
        event = self.source[idx]
        return self.get_traces_for_event(event)

    def _get_source(self, source):
        """Return the table containing the events.

        :param source: the *name* of the table.  If None, this method will
            try to find the original events table, even if the events were
            previously processed.
        :returns: table object

        """
        if source is None:
            if '_events' in self.group:
                source = self.group._events
            else:
                source = self.group.events
        else:
            source = self.data.get_node(self.group, source)
        return source

    def _check_destination(self, destination, overwrite):
        """Check if the destination is valid"""

        if destination == '_events':
            raise RuntimeError("The _events table is reserved for internal "
                               "use.  Choose another destination.")
        elif destination is None:
            destination = 'events'

        # If destination == source, source will be moved out of the way.  Don't
        # worry.  Otherwise, destination may not exist or will be overwritten
        if self.source.name != destination:
            if destination in self.group and not overwrite:
                raise RuntimeError("I will not overwrite previous results "
                                   "(unless you specify overwrite=True)")

        self.destination = destination

    def _clean_events_table(self):
        """Clean the events table.

        Remove duplicate events and sort the table by ext_timestamp.

        """
        events = self.source
        events_tablename = self.source.name

        enumerated_timestamps = \
            list(enumerate(events.col('ext_timestamp')))
        enumerated_timestamps.sort(key=operator.itemgetter(1))

        unique_sorted_ids = \
            self._find_unique_row_ids(enumerated_timestamps)

        new_events = self._replace_table_with_selected_rows(events,
                                                            unique_sorted_ids)
        self.source = new_events
        self._normalize_event_ids(new_events)

    def _find_unique_row_ids(self, enumerated_timestamps):
        """Find the unique row_ids from enumerated timestamps."""

        prev_timestamp = 0
        unique_sorted_ids = []
        for row_id, timestamp in enumerated_timestamps:
            if timestamp != prev_timestamp:
                # event is unique, so add it
                unique_sorted_ids.append(row_id)
            prev_timestamp = timestamp

        return unique_sorted_ids

    def _replace_table_with_selected_rows(self, table, row_ids):
        """Replace events table with selected rows.

        :param table: original table to be replaced.
        :param row_ids: row ids of the selected rows which should go in
            the destination table.

        """
        tmptable = self.data.create_table(self.group, 't__events',
                                          description=table.description)
        selected_rows = table.read_coordinates(row_ids)
        tmptable.append(selected_rows)
        tmptable.flush()
        self.data.rename_node(tmptable, table.name, overwrite=True)
        return tmptable

    def _normalize_event_ids(self, events):
        """Normalize event ids.

        After sorting, the event ids no longer correspond to the row
        number.  This can complicate finding the row id of a particular
        event.  This method will replace the event_ids of all events by
        the row id.

        :param events: the events table to normalize.

        """
        row_ids = range(len(events))
        events.modify_column(column=row_ids, colname='event_id')

    def _create_results_table(self):
        """Create results table containing the events."""

        self._tmp_events = self._create_empty_results_table()
        self._copy_events_into_table()

    def _create_empty_results_table(self):
        """Create empty results table with correct length."""

        if self.limit:
            length = self.limit
        else:
            length = len(self.source)

        if '_t_events' in self.group:
            self.data.remove_node(self.group, '_t_events')
        table = self.data.create_table(self.group, '_t_events',
                                       self.processed_events_description,
                                       expectedrows=length)

        for x in xrange(length):
            table.row.append()
        table.flush()

        return table

    def _copy_events_into_table(self):
        table = self._tmp_events
        source = self.source

        progressbar = self._create_progressbar_from_iterable(source.colnames)

        for col in progressbar(source.colnames):
            getattr(table.cols, col)[:self.limit] = getattr(source.cols,
                                                            col)[:self.limit]
        table.flush()

    def _store_results_from_traces(self):
        table = self._tmp_events

        timings = self.process_traces()

        # Assign values to full table, column-wise.
        for idx in range(4):
            col = 't%d' % (idx + 1)
            getattr(table.cols, col)[:] = timings[:, idx]
        table.flush()

    def process_traces(self, limit=None):
        """Process traces to yield pulse timing information."""

        if limit:
            self.limit = limit

        if self.limit:
            events = self.source.iterrows(stop=self.limit)
        else:
            events = self.source

        timings = self._process_traces_from_event_list(events,
                                                       length=self.limit)
        return timings

    def _process_traces_from_event_list(self, events, length=None):
        """Process traces from a list of events.

        This is the method looping over all events.

        :param events: an iterable of the events
        :param length: an indication of the number of events, for use as a
            progress bar.  Optional.

        """
        progressbar = self._create_progressbar_from_iterable(events, length)

        result = []
        for event in progressbar(events):
            timings = self._reconstruct_time_from_traces(event)
            result.append(timings)
        timings = np.array(result)

        return timings

    def _create_progressbar_from_iterable(self, iterable, length=None):
        """Create a progressbar object from any iterable."""

        if length is None:
            try:
                length = len(iterable)
            except TypeError:
                pass

        if length:
            return pb.ProgressBar(maxval=length, widgets=[pb.Percentage(),
                                                          pb.Bar(), pb.ETA()])
        else:
            # Cannot create progressbar, return no-op
            return lambda x: x

    def _reconstruct_time_from_traces(self, event):
        """Reconstruct arrival times for a single event.

        This method loops over the traces.

        :param event: row from the events table.
        :returns: arrival times in the detectors relative to trace start
                  in ns.

        """
        timings = []
        for baseline, pulseheight, trace_idx in zip(event['baseline'],
                                                    event['pulseheights'],
                                                    event['traces']):
            if pulseheight < 0:
                # retain -1, -999 status flags in timing
                timings.append(pulseheight)
            elif pulseheight < ADC_THRESHOLD:
                timings.append(-999)
            else:
                trace = self._get_trace(trace_idx)
                timings.append(self._reconstruct_time_from_trace(trace,
                                                                 baseline))
        timings = [time * ADC_TIME_PER_SAMPLE * 1e9
                   if time not in [-1, -999] else time
                   for time in timings]
        return timings

    def _get_trace(self, idx):
        """Returns a trace given an index into the blobs array.

        Decompress a trace from the blobs array.

        :param idx: index into the blobs array
        :returns: array of pulseheight values

        """
        blobs = self._get_blobs()

        trace = zlib.decompress(blobs[idx]).split(',')
        if trace[-1] == '':
            del trace[-1]
        trace = np.array([int(x) for x in trace])
        return trace

    def _get_blobs(self):
        return self.group.blobs

    def _reconstruct_time_from_trace(self, trace, baseline):
        """Reconstruct time of measurement from a trace.

        This method is doing the hard work.

        :param trace: array containing pulseheight values.
        :param baseline: baseline of the trace.
        :returns: index in trace for arrival time of first particle.

        """
        threshold = baseline + ADC_THRESHOLD
        value = self._first_above_threshold(trace, threshold)

        return value

    def _first_above_threshold(self, trace, threshold):
        """Find the first element in the list equal or above threshold

        If no element matches the condition -999 will be returned.

        """
        return next((i for i, x in enumerate(trace) if x >= threshold), -999)

    def _store_number_of_particles(self):
        """Store number of particles in the detectors.

        Process all pulseheights from the events and estimate the number
        of particles in each detector.

        """
        table = self._tmp_events

        n_particles = self._process_pulseintegrals()
        for idx in range(4):
            col = 'n%d' % (idx + 1)
            getattr(table.cols, col)[:] = n_particles[:, idx]
        table.flush()

    def _process_pulseintegrals(self):
        n_particles = []

        integrals = self.source.col('integrals')
        all_mpv = []
        for detector_integrals in integrals.T:
            if (detector_integrals < 0).all():
                all_mpv.append(np.nan)
            else:
                n, bins = np.histogram(detector_integrals,
                                       bins=np.linspace(0, 50000, 201))
                find_mpv = FindMostProbableValueInSpectrum(n, bins)
                mpv, is_fitted = find_mpv.find_mpv()
                if is_fitted:
                    all_mpv.append(mpv)
                else:
                    all_mpv.append(np.nan)
        all_mpv = np.array(all_mpv)

        for event in self.source[:self.limit]:
            pulseintegrals = event['integrals']
            # retain -1, -999 status flags
            pulseintegrals = np.where(pulseintegrals >= 0,
                                      pulseintegrals / all_mpv,
                                      pulseintegrals)
            # if mpv fit failed, value is nan.  Make it -999
            pulseintegrals = np.where(np.isnan(pulseintegrals), -999,
                                      pulseintegrals)
            n_particles.append(pulseintegrals)

        return np.array(n_particles)

    def _move_results_table_into_destination(self):
        if self.source.name == 'events':
            self.source.rename('_events')
            self.source = self.group._events

        if self.destination in self.group:
            self.data.remove_node(self.group, self.destination)
        self._tmp_events.rename(self.destination)

    def determine_detector_timing_offsets(self, timings_table='events'):
        """Determine the offsets between the station detectors."""

        table = self.data.get_node(self.group, timings_table)
        t2 = table.col('t2')

        gauss = lambda x, N, m, s: N * norm.pdf(x, m, s)
        bins = np.arange(-100 + 1.25, 100, 2.5)

        print "Determining offsets based on # events:",
        offsets = []
        for timings in 't1', 't3', 't4':
            timings = table.col(timings)
            dt = (timings - t2).compress((t2 >= 0) & (timings >= 0))
            print len(dt),
            y, bins = np.histogram(dt, bins=bins)
            x = (bins[:-1] + bins[1:]) / 2
            popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0., 10.))
            offsets.append(popt[1])
        print

        return [offsets[0]] + [0.] + offsets[1:]


class ProcessIndexedEvents(ProcessEvents):

    """Process a subset of events using an index.

    This is a subclass of :class:`ProcessEvents`.  Using an index, this
    class will only process a subset of events, thus saving time.  For
    example, this class can only process events making up a coincidence.

    """

    def __init__(self, data, group, indexes, source=None):
        """Initialize the class.

        :param data: the PyTables datafile
        :param group: the group containing the station data.  In normal
            cases, this is simply the group containing the events table.
        :param indexes: a list of indexes into the events table.
        :param source: the name of the events table.  Default: None,
            meaning the default name 'events'.

        """
        super(ProcessIndexedEvents, self).__init__(data, group, source)
        self.indexes = indexes

    def _store_results_from_traces(self):
        table = self._tmp_events

        timings = self.process_traces()

        for event, (t1, t2, t3, t4) in izip(table.itersequence(self.indexes),
                                            timings):
            event['t1'] = t1
            event['t2'] = t2
            event['t3'] = t3
            event['t4'] = t4
            event.update()

        table.flush()

    def process_traces(self):
        """Process traces to yield pulse timing information.

        This method makes use of the indexes to build a list of events.

        """
        events = self.source.itersequence(self.indexes)
        length = len(self.indexes)
        timings = self._process_traces_from_event_list(events, length=length)

        return timings

    def get_traces_for_indexed_event_index(self, idx):
        idx = self.indexes[idx]
        return self.get_traces_for_event_index(idx)


class ProcessEventsWithLINT(ProcessEvents):

    """Process events using LInear INTerpolation for arrival times.

    This is a subclass of :class:`ProcessEvents`.  Use a linear
    interpolation method to determine the arrival times of particles.

    """

    def _reconstruct_time_from_trace(self, trace, baseline):
        """Reconstruct time of measurement from a trace (LINT timings).

        This method is doing the hard work.

        :param trace: array containing pulseheight values.
        :param baseline: baseline of the trace
        :returns: arrival time

        """
        threshold = baseline + ADC_THRESHOLD
        i = self._first_above_threshold(trace, threshold)

        if i == 0:
            value = i
        elif not i == -999:
            x0, x1 = i - 1, i
            y0, y1 = trace[x0], trace[x1]
            value = 1. * (threshold - y0) / (y1 - y0) + x0
        else:
            value = -999

        return value


class ProcessIndexedEventsWithLINT(ProcessIndexedEvents,
                                   ProcessEventsWithLINT):

    """Process a subset of events using LInear INTerpolation.

    This is a subclass of :class:`ProcessIndexedEvents` and
    :class:`ProcessEventsWithLint`.

    """

    pass


class ProcessEventsWithoutTraces(ProcessEvents):

    """Process events without traces

    This is a subclass of :class:`ProcessEvents`.  Processing events
    without considering traces will invalidate the arrival time
    information.  However, for some analyses it is not necessary to obtain
    this information.  Ignoring the traces will then greatly decrease
    processing time and data size.

    """

    def _store_results_from_traces(self):
        """Fake storing results from traces."""

        pass


class ProcessIndexedEventsWithoutTraces(ProcessEventsWithoutTraces,
                                        ProcessIndexedEvents):

    """Process a subset of events without traces

    This is a subclass of :class:`ProcessIndexedEvents` and
    :class:`ProcessEventsWithoutTraces`.  Processing events without
    considering traces will invalidate the arrival time information.
    However, for some analyses it is not necessary to obtain this
    information.  Ignoring the traces will then greatly decrease
    processing time and data size.

    """

    pass


class ProcessEventsWithTriggerOffset(ProcessEvents):

    """Process events and reconstruct trigger time from traces

    The trigger times are stored in the columnt_trigger, they are
    relative to the start of traces, just like the t# columns.

    If no trigger can be found, possibly due to the data filter,
    a value of -999 will be entered.

    """

    def _store_results_from_traces(self):
        table = self._tmp_events

        timings = self.process_traces()

        # Assign values to full table, column-wise.
        for idx in range(4):
            col = 't%d' % (idx + 1)
            getattr(table.cols, col)[:] = timings[:, idx]
        getattr(table.cols, 't_trigger')[:] = timings[:, 4]
        table.flush()

    def _reconstruct_time_from_traces(self, event):
        """Reconstruct arrival times for a single event.

        This method loops over the traces.

        :param event: row from the events table.
        :returns: arrival times in the detectors and trigger time
                  relative to start of trace in ns

        """
        timings = []
        traces = []
        n_detectors = 4
        for baseline, pulseheight, trace_idx in zip(event['baseline'],
                                                    event['pulseheights'],
                                                    event['traces']):
            if pulseheight < 0:
                # retain -1, -999 status flags in timing
                timings.append(pulseheight)
                n_detectors -= 1
            elif pulseheight < ADC_THRESHOLD:
                timings.append(-999)
            else:
                trace = self._get_trace(trace_idx)
                traces.append(trace)
                timings.append(self._reconstruct_time_from_trace(trace,
                                                                 baseline))
        timings.append(self._reconstruct_trigger_time_from_traces(traces,
                                                                  n_detectors))
        timings = [time * ADC_TIME_PER_SAMPLE * 1e9
                   if time not in [-1, -999] else time
                   for time in timings]
        return timings

    def _reconstruct_trigger_time_from_traces(self, traces, n_detectors=None):
        """Reconstruct the moment of trigger from the traces

        The timestamp for an event is based on the moment of the
        trigger. This moment is reconstructed in this function. Using
        this result the arrival times in each detector can be corrected
        to be relative to the timestamp.

        This function assumes traces with no data filter applied and the
        default settings for triggers. The data filter in the HiSPARC
        DAQ can distort the peaks/pulses in traces, possibly preventing
        reconstruction.

        :param traces: the traces for an event.
        :param n_detectors: number of detectors, for trigger conditions.
        :returns: index in the trace for the trigger.

        """
        if not n_detectors:
            n_detectors = len(traces)
        if n_detectors not in [2, 4]:
            raise LookupError('Unsupported number of detectors')

        low_idx = []
        high_idx = []
        for trace in traces:
            low_idx.append(self._first_above_threshold(trace,
                                                       ADC_LOW_THRESHOLD))
            if n_detectors == 4 and not low_idx[-1] == -999:
                high_idx.append(self._first_above_threshold(
                    trace[low_idx[-1]:], ADC_HIGH_THRESHOLD))
                if not high_idx[-1] == -999:
                    high_idx[-1] += low_idx[-1]
        low_idx = [idx for idx in low_idx if not idx == -999]
        high_idx = [idx for idx in high_idx if not idx == -999]
        low_idx.sort()
        high_idx.sort()

        if n_detectors == 2 and len(low_idx) > 1:
            # Two low
            trigger_idx = low_idx[1]
        elif n_detectors == 4 and (len(low_idx) >= 3 or len(high_idx) >= 2):
            # Two high or three low
            if len(low_idx) < 3 and len(high_idx) >= 2:
                trigger_idx = high_idx[1]
            elif len(low_idx) >= 3 and len(high_idx) < 2:
                trigger_idx = low_idx[2]
            elif len(low_idx) >= 3 and len(high_idx) >= 2:
                trigger_idx = min([low_idx[2], high_idx[1]])
        else:
            # print "Trigger to low or filtered trace? max signals: ",
            # print [max(trace) for trace in traces]
            trigger_idx = -999

        return trigger_idx


class ProcessEventsFromSource(ProcessEvents):

    """Process HiSPARC events from a different source.

    This class is a subclass of ProcessEvents.  The difference is that in
    this class, the source and destination are assumed to be different
    files.  This also means that the source is untouched (no renaming of
    original event tables) and the destination is assumed to be empty.

    """

    def __init__(self, source_file, dest_file, source_group, dest_group):
        """Initialize the class.

        :param source_file: the PyTables source file
        :param dest_file: the PyTables dest file
        :param group_path: the pathname of the source (and destination)
            group

        """
        self.source_file = source_file
        self.dest_file = dest_file

        self.source_group = self.source_file.get_node(source_group)
        self.dest_group = self.dest_file.get_node(dest_group)

        self.source = self._get_source()

    def _get_source(self):
        """Return the table containing the events.

        :return: table object

        """
        if '_events' in self.source_group:
            source = self.source_group._events
        else:
            source = self.source_group.events
        return source

    def _check_destination(self, destination, overwrite):
        """Override method, the destination is empty"""
        pass

    def _replace_table_with_selected_rows(self, table, row_ids):
        """Replace events table with selected rows.

        :param table: original table to be replaced.
        :param row_ids: row ids of the selected rows which should go in
            the destination table.

        """
        new_events = self.dest_file.create_table(self.dest_group, '_events',
                                                 description=table.description)
        selected_rows = table.read_coordinates(row_ids)
        new_events.append(selected_rows)
        new_events.flush()
        return new_events

    def _create_empty_results_table(self):
        """Create empty results table with correct length."""

        if self.limit:
            length = self.limit
        else:
            length = len(self.source)

        table = self.dest_file.create_table(self.dest_group, 'events',
                                            self.processed_events_description,
                                            expectedrows=length)

        for x in xrange(length):
            table.row.append()
        table.flush()

        return table

    def _move_results_table_into_destination(self):
        """Override, destination is temporary table"""
        self.destination = self._tmp_events

    def _get_blobs(self):
        """Return blobs node"""

        return self.source_group.blobs

    def _create_progressbar_from_iterable(self, iterable, length=None):
        """Override method, do not show a progressbar"""

        return lambda x: x


class ProcessEventsFromSourceWithTriggerOffset(ProcessEventsFromSource,
                                               ProcessEventsWithTriggerOffset):

    """Process events from a different source and find trigger.

    This is a subclass of :class:`ProcessEventsFromSource` and
    :class:`ProcessEventsWithTriggerOffset`.  Processing events and
    finding the trigger time in the traces. And storing the results in a
    different file that the source.

    """

    pass
