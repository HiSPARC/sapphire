import zlib
from itertools import izip

import tables
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit
import progressbar as pb

from sapphire.storage import ProcessedHisparcEvent


ADC_THRESHOLD = 20
ADC_TIME_PER_SAMPLE = 2.5e-9


class ProcessEvents(object):
    def __init__(self, data, group, source=None):
        self.data = data
        self.group = data.getNode(group)
        self.source = self._get_source(source)

    def process_and_store_results(self, destination=None, overwrite=False,
                                  limit=None):
        self.limit = limit

        self._check_destination(destination, overwrite)

        self._create_results_table()
        self._store_results_from_traces()
        self._store_number_of_particles()
        self._move_results_table_into_destination()

    def get_traces_for_event(self, event):
        traces = []
        for idx in event['traces']:
            if not idx < 0:
                traces.append(self._get_trace(idx))

        # Make traces follow NumPy conventions
        traces = np.array(traces).T
        return traces

    def get_traces_for_event_index(self, idx):
        event = self.source[idx]
        return self.get_traces_for_event(event)

    def _get_source(self, source):
        if source is None:
            if '_events' in self.group:
                source = self.group._events
            else:
                source = self.group.events
        else:
            source = self.data.getNode(self.group, source)
        return source

    def _check_destination(self, destination, overwrite):
        if destination == '_events':
            raise RuntimeError("The _events table is reserved for internal use.  Choose another destination.")
        elif destination is None:
            destination = 'events'

        # If destination == source, source will be moved out of the way.  Don't
        # worry.  Otherwise, destination may not exist or will be overwritten
        if self.source.name != destination:
            if destination in self.group and not overwrite:
                raise RuntimeError("I will not overwrite previous results (unless you specify overwrite=True)")

        self.destination = destination

    def _create_results_table(self):
        self._tmp_events = self._create_empty_results_table()
        self._copy_events_into_table()

    def _create_empty_results_table(self):
        if self.limit:
            length = self.limit
        else:
            length = len(self.source)

        if '_t_events' in self.group:
            self.data.removeNode(self.group, '_t_events')
        table = self.data.createTable(self.group, '_t_events',
                                      ProcessedHisparcEvent,
                                      expectedrows=length)

        for x in xrange(length):
            table.row.append()
        table.flush()

        return table

    def _copy_events_into_table(self):
        table = self._tmp_events
        source = self.source

        progressbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])

        for col in progressbar(source.colnames):
            getattr(table.cols, col)[:self.limit] = getattr(source.cols,
                                                            col)[:self.limit]
        table.flush()

    def _store_results_from_traces(self):
        table = self._tmp_events

        timings = self.process_traces()

        for idx in range(4):
            col = 't%d' % (idx + 1)
            getattr(table.cols, col)[:] = timings[:, idx]
        table.flush()

    def process_traces(self, limit=None):
        """Process traces to yield pulse timing information"""

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
        progressbar = self._create_progressbar_from_iterable(events, length)

        result = []
        for event in progressbar(events):
            timings = self._reconstruct_time_from_traces(event)
            result.append(timings)
        result = np.array(result)

        # Replace NaN with -999, get timings in ns
        timings = np.where(np.isnan(result), -999, 1e9 * result)
        return timings

    def _create_progressbar_from_iterable(self, iterable, length=None):
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
        timings = []
        for pulseheight, trace_idx in zip(event['pulseheights'],
                                          event['traces']):
            if pulseheight < ADC_THRESHOLD:
                timings.append(np.nan)
            else:
                trace = self._get_trace(trace_idx)
                timings.append(self._reconstruct_time_from_trace(trace))
        return timings

    def _get_trace(self, idx):
        blobs = self.group.blobs

        trace = zlib.decompress(blobs[idx]).split(',')
        if trace[-1] == '':
            del trace[-1]
        trace = np.array([int(x) for x in trace])
        return trace

    def _reconstruct_time_from_trace(self, trace):
        """Reconstruct time of measurement from a trace"""

        t = trace[:100]
        baseline = np.mean(t)

        trace = trace - baseline
        threshold = ADC_THRESHOLD

        value = np.nan
        for i, t in enumerate(trace):
            if t >= threshold:
                value = i
                break

        return value * ADC_TIME_PER_SAMPLE

    def _store_number_of_particles(self):
        table = self._tmp_events

        n_particles = self.process_pulseheights()
        for idx in range(4):
            col = 'n%d' % (idx + 1)
            getattr(table.cols, col)[:] = n_particles[:, idx]
        table.flush()

    def process_pulseheights(self):
        #FIXME Much too simplistic!  Need fits

        n_particles = []

        for event in self.source[:self.limit]:
            n_particles.append(event['pulseheights'] / 380.)

        return np.array(n_particles)

    def _move_results_table_into_destination(self):
        if self.source.name == 'events':
            self.source.rename('_events')
            self.source = self.group._events

        if self.destination in self.group:
            self.data.removeNode(self.group, self.destination)
        self._tmp_events.rename(self.destination)

    def determine_detector_timing_offsets(self, timings_table='events'):
        table = self.data.getNode(self.group, timings_table)
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
    def __init__(self, data, group, indexes, source=None):
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
        events = self.source.itersequence(self.indexes)

        timings = self._process_traces_from_event_list(events,
                                                       length=len(self.indexes))
        return timings

    def get_traces_for_indexed_event_index(self, idx):
        idx = self.indexes[idx]
        return self.get_traces_for_event_index(idx)


class ProcessEventsWithLINT(ProcessEvents):
    def _reconstruct_time_from_trace(self, trace):
        """Reconstruct time of measurement from a trace (LINT timings)"""

        t = trace[:100]
        baseline = np.mean(t)

        trace = trace - baseline
        threshold = ADC_THRESHOLD

        # FIXME: apparently, there are a few bugs here. I see, in my
        # cluster reconstruction analysis, timings like -inf and
        # -something. Guesses: sometimes y0 == y1, and sometimes y1 < y0.

        value = np.nan
        for i, t in enumerate(trace):
            if t >= threshold:
                x0, x1 = i - 1, i
                y0, y1 = trace[x0], trace[x1]
                value = 1. * (threshold - y0) / (y1 - y0) + x0
                break

        return value * ADC_TIME_PER_SAMPLE


class ProcessIndexedEventsWithLINT(ProcessIndexedEvents, ProcessEventsWithLINT):
    pass
