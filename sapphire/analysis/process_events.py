import zlib

import tables
import numpy as np
import progressbar as pb


ADC_THRESHOLD = 20
ADC_TIME_PER_SAMPLE = 2.5e-9


class ProcessEvents():
    def __init__(self, data, group):
        self.data = data
        self.group = data.getNode(group)

    def process_traces(self, limit=None):
        """Process traces to yield pulse timing information"""

        if limit:
            events = self.group.events.iterrows(stop=limit)
        else:
            events = self.group.events

        return self.process_traces_from_event_list(events, length=limit)

    def process_traces_from_index(self, indexes):
        """Process traces to yield pulse timing information"""

        events = self.group.events

        return self.process_traces_from_event_list(events.itersequence(indexes),
                                                   length=len(indexes))

    def process_traces_from_event_list(self, events, length=None):
        if length is None:
            try:
                length = len(events)
            except TypeError:
                pass
        if length:
            progressbar = pb.ProgressBar(maxval=length,
                                         widgets=[pb.Percentage(), pb.Bar(),
                                                  pb.ETA()])
        else:
            progressbar = lambda x: x

        result = []
        for event in progressbar(events):
            timings = self.reconstruct_time_from_traces(event)
            result.append(timings)
        return 1e9 * np.array(result)

    def reconstruct_time_from_traces(self, event):
        timings = []
        for pulseheight, trace_idx in zip(event['pulseheights'],
                                          event['traces']):
            if pulseheight < ADC_THRESHOLD:
                timings.append(np.nan)
            else:
                trace = self.get_trace(trace_idx)
                timings.append(self.reconstruct_time_from_trace(trace))
        return timings

    def get_trace(self, idx):
        blobs = self.group.blobs

        trace = zlib.decompress(blobs[idx]).split(',')
        if trace[-1] == '':
            del trace[-1]
        trace = np.array([int(x) for x in trace])
        return trace

    def reconstruct_time_from_trace(self, trace):
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
