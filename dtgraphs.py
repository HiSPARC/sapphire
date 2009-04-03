import MySQLdb
import datetime
import time
import os

def calc_dt(x, kx, timeshift = 0):
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
        else:
            dt.append(v)
        ki += 1

    hist(dt, bins=100, range=(-1, 1), histtype='step',
         label="Shift %+g s" % timeshift)

    return dt

def read_from_db(limit):
    db = MySQLdb.connect('127.0.0.1', 'analysis', 'Data4analysis!', 'eventwarehouse', 3307)
    cursor = db.cursor()
    sql = "select date, time, nanoseconds, doublevalue from event join calculateddata using(event_id) join calculateddatatype using (calculateddatatype_id) where uploadcode='PH1' and station_id=601 limit %d" % limit
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
    kx = []; ky = []
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
        if t > ht[-1]:
            break
    return kx, ky

def do_timeshifts(shifts):
    for shift in shifts:
        print "Calculating dt's for timeshift", shift
        dt = calc_dt(x, kx, shift)

def finish_graph():
    legend()
    xlabel("Time difference (s)")
    ylabel("Counts")
    title("Nearest neighbour events for HiSPARC / KASCADE")

if __name__ == '__main__':
    os.environ['TZ'] = 'UTC'
    time.tzset()

    #result = read_from_db(100000)
    #d, x = calc_hisparc_times(result)
    #kx, ky = calc_kascade_times(d['time'])

    do_timeshifts([0, 1, -1, 13, -13, 13.2])
    finish_graph()
