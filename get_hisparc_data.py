""" Fetch events and eventdata from the eventwarehouse

    This module fetches data from the eventwarehouse, based on supplied
    start and end dates, limits and offsets. This data can then be passed
    on to the process_hisparc_data module.

"""
import MySQLdb
import datetime

def get_hisparc_data(station_id=601, start=None, stop=None, limit=None,
                     offset=None):
    """ Get events and eventdata from the eventwarehouse

    This function fetches events and eventdata from the eventwarehouse
    based on start and end dates, limits and offsets. You can pass this
    data on to the process_events function in the process_hisparc_data
    module.

    Arguments:
    station_id          the HiSPARC station id
    start               a datetime instance defining the start of the
                        search interval (inclusive)
    stop                a datetime instance defining the end of the search
                        interval (inclusive)
    limit               the maximum number of events
    offset              an offset in the total event list from which point
                        on a limit number of events is being selected

    """
    db = MySQLdb.connect('127.0.0.1', 'analysis', 'Data4analysis!',
                         'eventwarehouse', 3307)

    # get events from the event table
    events = get_hisparc_events(db, station_id, start, stop, limit, offset)

    if not events:
        return None, None

    # determine 'start' from the event data
    date = events[0][1]
    timedelta = events[0][2]
    start = datetime.datetime.combine(date, datetime.time()) + timedelta

    # determine 'stop' from the event data
    date = events[-1][1]
    timedelta = events[-1][2]
    stop = datetime.datetime.combine(date, datetime.time()) + timedelta

    # get the eventdata, where we don't select on event_ids, but rather
    # rely on 'start' and 'stop' instead.
    print "Time window: ", start, stop
    eventdata = get_hisparc_eventdata(db, station_id, start, stop)

    db.close()

    return events, eventdata

def get_hisparc_events(db, station_id, start=None, stop=None, limit=None,
                       offset=None):
    """ Get data from the eventd table

    This function fetches data from the event table, selected by using
    the start and stop datetime instances and possibly by using a limit and
    offset. This routine is generally not called directly. Instead, use
    get_hisparc_data which calls this routine and then continues to fetch
    the corresponding eventdata.

    Arguments:
    db                  a MySQLdb connection instance to the eventwarehouse
    station_id          the HiSPARC station id
    start               a datetime instance defining the start of the
                        search interval (inclusive)
    stop                a datetime instance defining the end of the search
                        interval (inclusive)
    limit               the maximum number of events
    offset              an offset in the total event list from which point on a
                        limit number of events is being selected

    """
    cursor = db.cursor()

    sql = "SELECT event_id, date, time, nanoseconds FROM event " \
          "WHERE station_id=%d AND eventtype_id=1 " % station_id
    if start:
        sql += "AND (date > '%s' OR (date = '%s' AND time >= '%s')) " \
               % (start.date(), start.date(),
                  start.time().strftime('%H:%M:%S'))
    if stop:
        sql += "AND (date < '%s' OR (date = '%s' AND time <= '%s')) " \
               % (stop.date(), stop.date(),
                  stop.time().strftime('%H:%M:%S'))
    sql += "ORDER BY date, time, nanoseconds "
    if limit:
        if offset:
            sql += "LIMIT %d,%d" % (offset, limit)
        else:
            sql += "LIMIT %d" % limit
    cursor.execute(sql)
    results = cursor.fetchall()
    return results    

def get_hisparc_eventdata(db, station_id, start=None, stop=None):
    """ Get data from the eventdata table

    This function fetches data from the eventdata table, selected by using
    the start and stop datetime instances. It doesn't select using a set of
    event_id's.  This is because we may have millions of events to select
    on. This inexact match will only be inexact at the end of the data,
    where we possibly include eventdata from extra events. The reason is
    that 'start' really matches the 'start' of the events, whereas 'stop'
    has to deal with possible 'limits' and makes sure to include all
    events, at the cost of possibly including more.

    It is important to realize that YOU are responsible for making this
    work by defining correct start and stop values.
    
    Of course, this is already done for you by the get_hisparc_data
    routine, which you should probably use.

    Arguments:
    db                  a MySQLdb connection instance to the eventwarehouse
    station_id          the HiSPARC station id
    start               a datetime instance defining the start of the
                        search interval (inclusive)
    stop                a datetime instance defining the end of the search
                        interval (inclusive)

    """
    cursor = db.cursor()

    sql  = "SELECT event_id, uploadcode, doublevalue " \
           "FROM event e JOIN calculateddata USING(event_id) " \
           "JOIN calculateddatatype USING(calculateddatatype_id) " \
           "WHERE station_id=%d AND e.eventtype_id=1 " % station_id
    sql += "AND uploadcode " \
           "IN ('PH1','PH2','PH3','PH4','IN1','IN2','IN3','IN4') "
    if start:
        sql += "AND (date > '%s' OR (date = '%s' AND time >= '%s')) " \
               % (start.date(), start.date(),
                  start.time().strftime('%H:%M:%S'))
    if stop:
        sql += "AND (date < '%s' OR (date = '%s' AND time <= '%s')) " \
               % (stop.date(), stop.date(),
                  stop.time().strftime('%H:%M:%S'))
    sql += "ORDER BY date, time, nanoseconds"
    cursor.execute(sql)
    results = cursor.fetchall()
    return results
