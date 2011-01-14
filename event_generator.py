import csv
import datetime
import struct

from numpy.random import randint, random
import MySQLdb


T0 = 1234567890 # 14 Feb 2009 00:31:30
H_SHIFT = 13.18 # HiSPARC timeshift

K_FILE = "generator-kascade.dat"

def generate_events(timespan, rate, reconstructed_fraction):
    """Generate HiSPARC and KASCADE synchronized events"""

    N = timespan * rate
    timestamps = randint(T0, T0 + timespan, N)
    nanoseconds = randint(0, int(1e9), N)

    with open(K_FILE, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
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
    t = datetime.datetime.fromtimestamp(ts)
    trace = 'xxx'
    l = len(trace)

    # num_devices = 2, use gps timestamp, ignore rest
    header = struct.pack('>2BBfBH5BH3L', 0, 0, 0, 0, 2, 2 * l, t.second,
                         t.minute, t.hour, t.day, t.month, t.year,
                         int(ns), 0, 0)
    master = struct.pack('>8H2L%ds%ds' % (l, l), 0, 0, 0, 0, 0, 0, 0, 0,
                         0, 0, trace, trace)
    slave = struct.pack('>8H2L%ds%ds' % (l, l), 0, 0, 0, 0, 0, 0, 0, 0, 0,
                        0, trace, trace)

    msg = header + master + slave

    cursor.execute("INSERT INTO message (device_id, message) VALUES "
                   "(601, %s)", (msg,))

def store_kascade_event(writer, ts, ns):
    ns = ns - (ns % 200)
    writer.writerow((0, 0, ts, ns, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0))


if __name__ == '__main__':
    generate_events(3600, 4., .1)
