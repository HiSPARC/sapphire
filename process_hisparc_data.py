import datetime
import time
import os

def process_events(events, table):
    tablerow = table.row

    print "Calculating and loading event timestamps..."
    for row in events:
        event_id = row[0]
        date = row[1]
        timedelta = row[2]
        nanoseconds = row[3]

        t = datetime.datetime.combine(date, datetime.time()) + timedelta
        timestamp = time.mktime(t.utctimetuple())

        tablerow['event_id'] = event_id
        tablerow['timestamp'] = timestamp
        tablerow['nanoseconds'] = nanoseconds
        tablerow.append()

    table.flush()

def process_eventdata(eventdata, table):
    print "Loading event data..."
    table_index = table.col('event_id').tolist()
    for row in eventdata:
        event_id = row[0]
        uploadcode = row[1]
        value = row[2]

        if uploadcode[:2] == 'PH':
            key = 'pulseheights'
        elif uploadcode[:2] == 'IN':
            key = 'integrals'
        else:
            continue

        try:
            tidx = table_index.index(event_id)
        except ValueError:
            continue
        else:
            data = table.read(tidx, field=key)[0]
            idx = int(uploadcode[2]) - 1 
            data[idx] = value
            table.modifyColumn(tidx, column=data, colname=key)

    table.flush()
        
if __name__ == '__main__':
    os.environ['TZ'] = 'UTC'
    time.tzset()
