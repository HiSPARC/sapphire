""" Read data from a raw data file and save timestamps as CSV

    This example script reads event data from a raw data file, extracts
    the timestamps, reformats them in date, time strings and saves them,
    along with nanoseconds, into a CSV file.

    This script was used to provide a high school student with detector
    data for his analysis.
"""

import tables
import datetime
import csv
import zlib

def get_traces(groupnode, traces_idx):
    """Get traces from groupnode and reference"""

    traces_array = groupnode.traces
    traces = []
    for idx in traces_idx:
        trace = zlib.decompress(traces_array[idx])
        trace = [int(x) * -.57 + 114.0 for x in trace.split(',')[:-1]]
        trace = ','.join(['%d' % x for x in trace])
        traces.append(trace)

    return traces


if __name__ == '__main__':
    # Read data from file
    with tables.openFile('ettyhillesum.h5', 'r') as datafile:
        data = [(x['timestamp'], x['nanoseconds'], x['pulseheights'],
                 x['integrals'], x['n_peaks'], x['traces']) for x in
                datafile.root.hisparc.ettyhillesum.events]

        # Reformat data
        datalist = []
        for timestamp, nanoseconds, pulseheights, integrals, n_peaks, \
            traces in data:
            # Format timestamps into date and time strings
            d_t = datetime.datetime.fromtimestamp(timestamp)
            date = d_t.date().isoformat()
            time = d_t.time().isoformat()
            # Convert pulseheights / integrals into mV units (only 2 master
            # channels)
            pulseheights = [x * .57 for x in pulseheights[:2]]
            integrals = [x * .57 for x in integrals[:2]]
            ph1, ph2 = pulseheights
            in1, in2 = integrals
            np1, np2 = n_peaks[:2]
            tr1, tr2 = get_traces(datafile.root.hisparc.ettyhillesum,
                                  traces[:2])
            # Save data values
            datalist.append((date, time, nanoseconds, ph1, ph2, in1, in2,
                             np1, np2, tr1, tr2))

    # Write data into a CSV file
    with open('ettyhillesum.csv', 'w') as file:
        writer = csv.writer(file, dialect='excel-tab')
        writer.writerows(datalist)
