""" Process HiSPARC traces

    This module can be used analyse (raw) traces.

"""
from numpy import around, mean, sign


FILTER_THRESHOLD = 20


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
