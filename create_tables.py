import tables

class HisparcEvent(tables.IsDescription):
    event_id = tables.Int64Col(pos=0)
    timestamp =  tables.Time32Col(pos=1)
    nanoseconds = tables.UInt32Col(pos=2)
    pulseheights = tables.Int16Col(shape=4, dflt=-9999, pos=3)
    integrals = tables.Int32Col(shape=4, dflt=-9999, pos=4)

data = tables.openFile('data.h5', 'w', 'HiSPARC / KASCADE data')
hisparc = data.createGroup('/', 'hisparc', 'HiSPARC data')
table = data.createTable(hisparc, 'events', HisparcEvent, 'HiSPARC events')
data.close()
