""" Process HiSPARC traces

    This module can be used analyse (raw) traces. It implements the same
    algorithms as are implemented in the HiSPARC DAQ.

    The :class:`MeanFilter` is meant to mimic the filter in the HiSPARC DAQ.
    It is reproduced here to make it easy to read the algorithm.

"""
from numpy import around, convolve, ones, where
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

    """Reconstruct trace observables

    If one wants to reconstruct trace observables from existing data some
    caveats apply. If the station applied the Mean Filter the trace values
    will no longer match the raw values used to determine the observables
    on the station. Additionally, if data reduction was active the trace
    may be missing samples without a significant signal, this complicates
    the determination of the baseline. Moreover, data reduction uses the
    BASELINE_THRESHOLD to determine what signals to keep, so tiny pulses
    may be removed making it impossible to reconstruct tiny pulseheights.

    The (default) value of BASELINE_THRESHOLD is different for the
    HiSPARC DAQ prior to v4 and also for PySPARC. Those use 25 ADC as
    threshold.

    """

    def __init__(self, traces):
        self.traces = traces
        self.n = self.traces.shape[1]
        self.slave = []
        if self.n == 2:
            self.slave = [-1, -1]
        elif not self.n == 4:
            raise Exception('Unsupported number of detectors')

    @lazy
    def baselines(self):
        """Mean value of the first 50 samples of the trace

        Usually this value is either around 200 or 30.

        :return: the baseline in ADC count.

        """
        baselines = around(self.traces[:50].mean(axis=0)).astype('int')
        return baselines.tolist() + self.slave

    @lazy
    def std_dev(self):
        """Standard deviation of the first 50 samples of the trace

        :return: the standard deviation in milli ADC count.

        """
        std_dev = around(self.traces[:50].std(axis=0) * 1000).astype('int')
        return std_dev.tolist() + self.slave

    @lazy
    def pulseheights(self):
        """Maximum peak to baseline value in trace

        :return: the pulseheights in ADC count.

        """
        pulseheights = self.traces.max(axis=0) - self.baselines[:self.n]
        return pulseheights.tolist() + self.slave

    @lazy
    def integrals(self):
        """Integral of trace for all values over threshold

        The threshold is defined by BASELINE_THRESHOLD

        :return: the pulse integral in ADC count * sample.

        """
        threshold = BASELINE_THRESHOLD
        integrals = where(self.traces - self.baselines[:self.n] > threshold,
                          self.traces - self.baselines[:self.n], 0).sum(axis=0)
        return integrals.tolist() + self.slave


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

        moving_average = convolve(trace, ones(4) / 4)

        filtered_trace = []
        local_mean = moving_average[3]

        if all([abs(v - local_mean) <= self.threshold for v in trace[:4]]):
            filtered_trace.extend([int(around(local_mean))] * 4)
        else:
            filtered_trace.extend(trace[:4])

        for i in xrange(4, len(trace)):
            local_mean = moving_average[i]
            if abs(trace[i] - trace[i - 1]) > 2 * self.threshold:
                filtered_trace.append(trace[i])
            elif (trace[i] > local_mean) == (trace[i - 1] > local_mean):
                # Both values on same side of the local_mean
                filtered_trace.append(trace[i])
            elif abs(trace[i] - local_mean) > self.threshold:
                filtered_trace.append(trace[i])
            else:
                filtered_trace.append(int(around(local_mean)))

        return filtered_trace

    def mean_filter_without_threshold(self, trace):
        """The mean filter in case use_threshold is False"""

        moving_average = convolve(trace, ones(4) / 4)
        local_mean = moving_average[3]
        filtered_trace = [int(around(local_mean))] * 4

        for i in xrange(4, len(trace)):
            local_mean = moving_average[i]
            if (trace[i] > local_mean) == (trace[i - 1] > local_mean):
                # Both values on same side of the local_mean
                filtered_trace.append(trace[i])
            else:
                filtered_trace.append(int(around(local_mean)))

        return filtered_trace
