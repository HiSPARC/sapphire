"""Check particle arriving time code.

Arne found a problem with the particle arrival time code
(:meth:`sapphire.analysis.process_events.ProcessEvents._reconstruct_time_from_trace`
and
:meth:`sapphire.analysis.process_events.ProcessEventsWithLINT._reconstruct_time_from_trace`).
The code was only valid for full traces (with 'reduce data' turned off)
because 'reduce data' did not provide a leader of 100 baseline samples.
The code is fixed now and this test makes sure you're actually using the
new code.

"""
import numpy as np
import tables

from sapphire.analysis import process_events


def main():
    # Create bogus data file
    data = tables.open_file('test.h5', 'w')
    data.create_array('/', 'events', np.arange(10))

    # Create bogus trace
    trace = np.arange(0, 150)
    trace *= 10
    trace += 5

    # Get timings from trace
    process = process_events.ProcessEvents(data, '/')
    t0 = process._reconstruct_time_from_trace(trace, 0.)
    process = process_events.ProcessEventsWithLINT(data, '/')
    t1 = process._reconstruct_time_from_trace(trace, 0.)
    t0 /= 2.5e-9
    t1 /= 2.5e-9
    print t0, t1

    # Assertions
    assert abs(t0 - 2.0) < 1e-6
    assert abs(t1 - 1.5) < 1e-6


if __name__ == '__main__':
    main()
