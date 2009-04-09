import MySQLdb
import datetime
import time
import os

def calc_dt(x, kx, timeshift = 0):
    matches = []

    kx = [u + timeshift for u in kx]
    dt = []
    
    for ki in range(len(kx)):
        if kx[ki] > x[0]:
            break

    i = 1
    while True:
        try:
            kv = kx[ki]
            hv = x[i]
            hpv = x[i-1]
        except IndexError:
            break
        if kv > hv:
            i += 1
            continue
        u = kv - hv
        v = kv - hpv
        if abs(u) < abs(v):
            dt.append(u)
            matches.append((u, i, ki))
        else:
            dt.append(v)
            matches.append((v, i-1, ki))
        ki += 1

    hist(dt, bins=100, range=(-1, 1), histtype='step',
         label="Shift %+g s" % timeshift)

    return matches

def read_from_db(limit):
    db = MySQLdb.connect('127.0.0.1', 'analysis', 'Data4analysis!',
                         'eventwarehouse', 3307)
    cursor = db.cursor()
    sql = "SELECT date, time, nanoseconds, doublevalue " \
          "FROM event " \
          "JOIN calculateddata USING(event_id) " \
          "JOIN calculateddatatype USING (calculateddatatype_id) " \
          "WHERE uploadcode='IN2' " \
          "AND station_id=601 " \
          "LIMIT %d" % limit
    cursor.execute(sql)
    result = cursor.fetchall()
    return result

def calc_hisparc_times(result):
    d = {}; d['time'] = []; d['nanoseconds'] = []; d['ph'] = []
    for r in result:
        t = datetime.datetime.combine(r[0], datetime.time()) + r[1]
        ts = time.mktime(t.utctimetuple())
        d['time'].append(ts)
        d['nanoseconds'].append(r[2]); d['ph'].append(r[3])
        
    x = []
    for i in range(len(d['time'])):
        x.append(d['time'][i] + d['nanoseconds'][i] / 1e9)

    return d, x

def calc_kascade_times(ht):
    kx = []; ky = []; kyy = []
    f = open('kascade-time.dat', 'r')
    while True:
        line = f.readline()
        data = line.split(' ')
        data = [float(v) for v in data]
        t = data[0] + data[1] / 1e9
        if t < ht[0]:
            continue
        kx.append(t)
        ky.append(data[4])
        kyy.append(data[3])
        if t > ht[-1]:
            break
    return kx, ky, kyy

def do_timeshifts(shifts):
    for shift in shifts:
        print "Calculating dt's for timeshift", shift
        matches = calc_dt(x, kx, shift)
    return matches

def finish_graph():
    legend()
    xlabel("Time difference (s)")
    ylabel("Counts")
    title("Nearest neighbour events for HiSPARC / KASCADE")

if __name__ == '__main__':
    os.environ['TZ'] = 'UTC'
    time.tzset()

    result = read_from_db(1000000)
    d, x = calc_hisparc_times(result)
    kx, ky, kyy = calc_kascade_times(d['time'])

    matches = do_timeshifts([13.180212926864623])
    finish_graph()

    matches = [v for v in matches if abs(v[0]) < 1e-5]
