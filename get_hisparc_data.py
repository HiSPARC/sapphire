import MySQLdb
import datetime

def get_hisparc_data(station_id=601, start=None, stop=None, limit=None,
                     offset=None):
    db = MySQLdb.connect('127.0.0.1', 'analysis', 'Data4analysis!',
                         'eventwarehouse', 3307)

    events = get_hisparc_events(db, station_id, start, stop, limit, offset)

    if not start:
        date = events[0][1]
        timedelta = events[0][2]
        start = datetime.datetime.combine(date, datetime.time()) + timedelta
    if not stop:
        date = events[-1][1]
        timedelta = events[-1][2]
        stop = datetime.datetime.combine(date, datetime.time()) + timedelta

    eventdata = get_hisparc_eventdata(db, station_id, start, stop)

    return events, eventdata

def get_hisparc_events(db, station_id, start, stop, limit, offset):
    cursor = db.cursor()

    sql = "SELECT event_id, date, time, nanoseconds FROM event " \
          "WHERE station_id=%d AND eventtype_id=1 " % station_id
    if start:
        sql += "AND date >= '%s' AND time >= '%s' " \
               % (start.date(), start.time().strftime('%H:%M:%S'))
    if stop:
        sql += "AND date <= '%s' AND time <= '%s' " \
               % (stop.date(), stop.time().strftime('%H:%M:%S'))
    sql += "ORDER BY date, time, nanoseconds "
    if limit:
        if offset:
            sql += "LIMIT %d,%d" % (offset, limit)
        else:
            sql += "LIMIT %d" % limit
    cursor.execute(sql)
    results = cursor.fetchall()
    return results    

def get_hisparc_eventdata(db, station_id, start, stop):
    cursor = db.cursor()

    sql  = "SELECT event_id, uploadcode, doublevalue " \
           "FROM event e JOIN calculateddata USING(event_id) " \
           "JOIN calculateddatatype USING(calculateddatatype_id) " \
           "WHERE station_id=%d AND e.eventtype_id=1 " % station_id
    sql += "AND date >= '%s' AND time >= '%s' " \
           % (start.date(), start.time().strftime('%H:%M:%S'))
    sql += "AND date <= '%s' AND time <= '%s' " \
           % (stop.date(), stop.time().strftime('%H:%M:%S'))
    sql += "AND uploadcode " \
           "IN ('PH1','PH2','PH3','PH4','IN1','IN2','IN3','IN4') " \
           "ORDER BY date, time, nanoseconds"
    cursor.execute(sql)
    results = cursor.fetchall()
    return results    
