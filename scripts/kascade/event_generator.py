import csv
import datetime
import struct

from numpy.random import randint, random
import MySQLdb

from sapphire.transformations import clock


T0 = 1234567890 # 14 Feb 2009 00:31:30
H_SHIFT = 13.18 # HiSPARC timeshift

K_FILE = "generator-kascade.dat"


def generate_events(timespan, rate, reconstructed_fraction):
    """Generate HiSPARC and KASCADE synchronized events"""

    N = timespan * rate
    timestamps = randint(T0, T0 + timespan, N)
    nanoseconds = randint(0, int(1e9), N)

    timestamps, nanoseconds = zip(*sorted(zip(timestamps, nanoseconds)))

    with open(K_FILE, 'w') as f:
        writer = csv.writer(f, delimiter=' ')
        db = MySQLdb.connect(user='buffer', passwd='Buffer4hisp!',
                             db='buffer')
        cursor = db.cursor()
        cursor.execute("DELETE FROM message")
        db.commit()

        for ts, ns in zip(timestamps, nanoseconds):
            store_hisparc_event(cursor, ts, ns)
            if random() < reconstructed_fraction:
                store_kascade_event(writer, ts, ns)
    db.commit()
    db.close()


def store_hisparc_event(cursor, ts, ns):
    t = datetime.datetime.utcfromtimestamp(clock.utc_to_gps(ts))
    trace = 'xxx'
    l = len(trace)

    msg = struct.pack(">BBHBBBBBHIiHhhhhhhii%ds%dshhhhhhii%ds%ds" %
                      (l, l, l, l),
                      2,        # central database
                      2,        # number of devices
                      l * 2,    # length of two traces
                      t.second,
                      t.minute,
                      t.hour,
                      t.day,
                      t.month,
                      t.year,
                      int(ns),
                      0,        # SLVtime
                      0,        # Trigger pattern
                      0,        # baseline1
                      0,        # baseline2
                      0,        # npeaks1
                      0,        # npeaks2
                      0,        # pulseheight1
                      0,        # pulseheight2
                      0,        # integral1
                      0,        # integral2
                      trace, trace,
                      0,        # baseline3
                      0,        # baseline4
                      0,        # npeaks3
                      0,        # npeaks4
                      0,        # pulseheight3
                      0,        # pulseheight4
                      0,        # integral3
                      0,        # integral4
                      trace, trace)

    cursor.execute("INSERT INTO message (device_id, message) VALUES "
                   "(601, %s)", (msg,))


def store_kascade_event(writer, ts, ns):
    ns = ns - (ns % 200)
    writer.writerow((0, 0, ts, ns, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0))


if __name__ == '__main__':
    generate_events(86400, 4., .1)
