""" Analyse HiSPARC traces

    Using this module you can analyse HiSPARC traces.

"""
import zlib
import numpy as np


def get_traces(traces_array, indexes):
    """Get traces from the traces table with given indexes

    This function will process compressed traces and transform them in
    an array of values in millivolts.

    :param traces_array: the PyTables traces_array which contains the
        compressed traces
    :param indexes: a list or array containing the indexes of the traces
        into the traces_array

    Example usage::

        >>> import tables
        >>> data = tables.openFile('test.h5', 'a')
        >>> traces_i = data.root.hisparc.station505.events[0]['traces']
        >>> from hisparc.analysis.traces import get_traces
        >>> traces = get_traces(data.root.hisparc.station505.traces,
        ... traces_i)

    Now you can plot the traces, for example.

    """
    traces = []

    for i in indexes:
        trace = zlib.decompress(traces_array[i]).split(',')

        # See if the last value of the trace is empty.  if so, remove it
        if trace[-1] == '':
            trace = trace[:-1]

        traces.append(np.array(map(adc_to_mv, trace)))

    return traces


def adc_to_mv(value):
    """Transform an ADC value into millivolts"""

    #FIXME
    return (int(value) - 200) * -.57
