import tables

class HisparcEvent(tables.IsDescription):
    event_id = tables.UInt64Col(pos=0)
    timestamp = tables.Time32Col(pos=1)
    nanoseconds = tables.UInt32Col(pos=2)
    pulseheights = tables.Int16Col(shape=4, dflt=-9999, pos=3)
    integrals = tables.Int32Col(shape=4, dflt=-9999, pos=4)

class KascadeEvent(tables.IsDescription):
    run_id = tables.IntCol(pos=0)
    event_id = tables.Int64Col(pos=1)
    timestamp = tables.Time32Col(pos=2)
    nanoseconds = tables.UInt32Col(pos=3)
    energy = tables.FloatCol(pos=4)
    core_pos = tables.FloatCol(shape=2, pos=5)
    zenith = tables.FloatCol(pos=6)
    azimuth = tables.FloatCol(pos=7)
    Num_e = tables.FloatCol(pos=8)
    Num_mu = tables.FloatCol(pos=9)
    dens_e = tables.FloatCol(shape=4, pos=10)
    dens_mu = tables.FloatCol(shape=4, pos=11)
    P200 = tables.FloatCol(pos=12)
    T200 = tables.FloatCol(pos=13)


if __name__ == '__main__':
    data = tables.openFile('data.h5', 'w', 'HiSPARC / KASCADE data')
    hisparc = data.createGroup('/', 'hisparc', 'HiSPARC data')
    kascade = data.createGroup('/', 'kascade', 'KASCADE data')
    data.createTable(hisparc, 'events', HisparcEvent, 'HiSPARC events')
    data.createTable(kascade, 'events', KascadeEvent, 'KASCADE events')
    data.close()
