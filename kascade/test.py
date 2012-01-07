from __future__ import division

import tables
import time

from sapphire.analysis import process_events


N = 1000


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'r')

    cls = process_events.ProcessEvents(
                data, '/hisparc/cluster_kascade/station_601')

    t0 = time.time()
    ts = cls.process_traces(limit=N)
    dt = time.time() - t0
    guess_dt = (len(data.root.hisparc.cluster_kascade.station_601.events) /
                N * dt)
    print "Processing %d events: %.2f s" % (N, dt)
    print "Estimate for total processing time: %d s" % guess_dt

    t0 = time.time()
    indexes = data.root.kascade.c_index[:,1][:N]
    ts = cls.process_traces_from_index(indexes)
    dt = time.time() - t0
    guess_dt = len(data.root.kascade.c_index) / N * dt
    print "Processing %d events: %.2f s" % (N, dt)
    print "Estimate for total processing time: %d s" % guess_dt
