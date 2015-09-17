""" Process HiSPARC traces

    This module can be used analyse (raw) traces. It implements the same
    algorithms as are implemented in the HiSPARC DAQ.

    The :class:`MeanFilter` is meant to mimic the filter in the HiSPARC DAQ.
    It is reproduced here to make it easy to read the algorithm.

"""
from numpy import around, mean, sign, std
from lazy import lazy


FILTER_THRESHOLD = 10  # default -6 mV
BASELINE_THRESHOLD = 17  # default -10 mV
ADC_TIME_PER_SAMPLE = 2.5  # in ns

# Trigger windows
PRE_TRIGGER = 400  # samples, i.e. 1000 ns
TRIGGER = 600  # samples, i.e. 1500 ns
POST_TRIGGER = 1400  # samples, i.e. 3500 ns

# Trigger thresholds
# HiSPARC II with baseline at 200 ADC counts.
ADC_LOW_THRESHOLD = 253
ADC_HIGH_THRESHOLD = 323

# HiSPARC III with baseline at 30 ADC counts (DAQ v4).
ADC_LOW_THRESHOLD_III = 82
ADC_HIGH_THRESHOLD_III = 150


class TraceObservables(object):

    """Reconstruct trace observables"""

    def __init__(self, traces):
        self.traces = traces

    @lazy
    def baselines(self):
        """Mean value of the first 100 samples of the trace"""

        return [int(around(mean(t[:100]))) for t in self.traces]

    @lazy
    def std_dev(self):
        """Standard deviation of the first 100 samples of the trace"""

        return [int(round(std(t[:100]))) for t in self.traces]

    @lazy
    def pulseheights(self):
        """Maximum peak to baseline value in trace"""

        return [max(t) - b for t, b in zip(self.traces, self.baselines)]

    @lazy
    def integrals(self):
        """Integral of trace for all values over threshold

        The threshold is defined by BASELINE_THRESHOLD

        """
        return [sum(v - b for v in trace if v - b > BASELINE_THRESHOLD)
                for trace, b in zip(self.traces, self.baselines)]


class TriggerReconstruction(object):

    """Reconstruct the sample in a trace where the trigger occurred"""

    pass


class MeanFilter(object):

    """Filter raw traces

    This class replicates the behavior of the Mean_Filter.vi in the HiSPARC
    DAQ. A sawtooth-like pattern may be visible in traces due to ADC
    misalignment and the synchronization pulse (?). This filter removes such
    small oscillations but keeps significant pulses.

    .. warning::

        This is a destructive algorithm which removes important information
        from the data.

    Verified, by eye, with the LabView VI and Bachelor thesis Oostenbrugge2014.

    """

    def __init__(self, use_threshold=True, threshold=FILTER_THRESHOLD):
        """Initialize the class.

        :param use_threshold: use a threshold when filtering traces.
        :param threshold: value of the threshold to use.

        """
        if use_threshold:
            self.filter = self.mean_filter_with_threshold
            self.threshold = threshold
        else:
            self.filter = self.mean_filter_without_threshold

    def filter_traces(self, raw_traces):
        """Apply the mean filter to multiple traces"""

        return [self.filter_trace(raw_trace) for raw_trace in raw_traces]

    def filter_trace(self, raw_trace):
        """Apply the mean filter to a single trace

        First separate the even and odd ADC traces, filter each separately.
        Then recombine them and pass the entire trace through the filter.

        """
        even_trace = raw_trace[::2]
        odd_trace = raw_trace[1::2]

        filtered_even = self.filter(even_trace)
        filtered_odd = self.filter(odd_trace)

        recombined_trace = [v
                            for eo in zip(filtered_even, filtered_odd)
                            for v in eo]
        filtered_trace = self.filter(recombined_trace)
        return filtered_trace

    def mean_filter_with_threshold(self, trace):
        """The mean filter in case use_threshold is True"""

        filtered_trace = []
        local_mean = mean(trace[:4])

        if all([abs(v - local_mean) <= self.threshold for v in trace[:4]]):
            filtered_trace.extend([int(around(local_mean))] * 4)
        else:
            filtered_trace.extend(trace[:4])

        for i in xrange(4, len(trace)):
            if abs(trace[i] - trace[i - 1]) > 2 * self.threshold:
                filtered_trace.append(trace[i])
            elif (sign(trace[i] - local_mean) ==
                  sign(trace[i - 1] - local_mean)):
                # Both values on same side of the local_mean
                filtered_trace.append(trace[i])
            elif abs(trace[i] - local_mean) > self.threshold:
                filtered_trace.append(trace[i])
            else:
                filtered_trace.append(int(around(local_mean)))

        return filtered_trace

    def mean_filter_without_threshold(self, trace):
        """The mean filter in case use_threshold is False"""

        local_mean = mean(trace[:4])
        filtered_trace = [int(around(local_mean))] * 4

        for i in xrange(4, len(trace)):
            local_mean = mean(trace[i - 4:i])
            if sign(trace[i] - local_mean) == sign(trace[i - 1] - local_mean):
                # Both values on same side of the local_mean
                filtered_trace.append(trace[i])
            else:
                filtered_trace.append(int(around(local_mean)))
        return filtered_trace
