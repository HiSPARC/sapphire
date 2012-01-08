import zlib

import tables
import numpy as np
import progressbar as pb

from sapphire.storage import ProcessedHisparcEvent


ADC_THRESHOLD = 20
ADC_TIME_PER_SAMPLE = 2.5e-9


class ProcessEvents(object):
    def __init__(self, data, group, limit=None, overwrite=False):
        self.data = data
        self.group = data.getNode(group)
        self.limit = limit
        self.overwrite = overwrite

    def process_and_store_results(self):
        if '_events' in self.group:
            if not self.overwrite:
                raise RuntimeError("I found an _events node.  Will not overwrite previous results")
            else:
                self.group.events.remove()
                self.group._events.rename('events'
                                          )
        self._create_results_table()
        self._store_results_from_traces()
        self._store_number_of_particles()
        self._move_results_table()

    def _create_results_table(self):
        self._tmp_events = self._create_empty_results_table()
        self._copy_events_into_table()

    def _create_empty_results_table(self):
        table = self.data.createTable(self.group, '_t_events',
                                      ProcessedHisparcEvent)
        if self.limit:
            length = self.limit
        else:
            length = len(self.group.events)

        for x in xrange(length):
            table.row.append()
        table.flush()

        return table

    def _copy_events_into_table(self):
        table = self._tmp_events
        events = self.group.events

        for col in events.colnames:
            getattr(table.cols, col)[:self.limit] = getattr(events.cols,
                                                            col)[:self.limit]
        table.flush()

    def _store_results_from_traces(self):
        table = self._tmp_events

        timings = self.process_traces()
        for idx in range(4):
            col = 't%d' % (idx + 1)
            getattr(table.cols, col)[:] = timings[:, idx]
        table.flush()

    def process_traces(self):
        """Process traces to yield pulse timing information"""

        if self.limit:
            events = self.group.events.iterrows(stop=self.limit)
        else:
            events = self.group.events

        timings = self._process_traces_from_event_list(events,
                                                       length=self.limit)
        return timings

    def _process_traces_from_event_list(self, events, length=None):
        progressbar = self._create_progressbar_from_iterable(events, length)

        result = []
        for event in progressbar(events):
            timings = self._reconstruct_time_from_traces(event)
            result.append(timings)
        return 1e9 * np.array(result)

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

        for event in self.group.events[:self.limit]:
            n_particles.append(event['pulseheights'] / 380.)

        return np.array(n_particles)

    def _move_results_table(self):
        self.group.events.rename('_events')
        self._tmp_events.rename('events')


class ProcessIndexedEvents(ProcessEvents):
    def __init__(self, data, group, indexes):
        super(ProcessIndexedEvents, self).__init__(data, group)
        self.indexes = indexes

    def process_traces(self):
        events = self.group.events.itersequence(self.indexes)

        timings = self._process_traces_from_event_list(events,
                                                       length=len(self.indexes))
        return timings
