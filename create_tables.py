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

class Coincidence(tables.IsDescription):
    hisparc_event_id = tables.UInt64Col(pos=0)
    kascade_event_id = tables.UInt64Col(pos=1)
    hisparc_timestamp = tables.Time32Col(pos=2)
    hisparc_nanoseconds = tables.UInt32Col(pos=3)
    hisparc_pulseheights = tables.Int16Col(shape=4, dflt=-9999, pos=4)
    hisparc_integrals = tables.Int32Col(shape=4, dflt=-9999, pos=5)
    kascade_timestamp = tables.Time32Col(pos=6)
    kascade_nanoseconds = tables.UInt32Col(pos=7)
    kascade_energy = tables.FloatCol(pos=8)
    kascade_core_pos = tables.FloatCol(shape=2, pos=9)
    kascade_zenith = tables.FloatCol(pos=10)
    kascade_azimuth = tables.FloatCol(pos=11)
    kascade_Num_e = tables.FloatCol(pos=12)
    kascade_Num_mu = tables.FloatCol(pos=13)
    kascade_dens_e = tables.FloatCol(shape=4, pos=14)
    kascade_dens_mu = tables.FloatCol(shape=4, pos=15)
    kascade_P200 = tables.FloatCol(pos=16)
    kascade_T200 = tables.FloatCol(pos=17)


if __name__ == '__main__':
    print "Creating a new PyTables data file... ",
    data = tables.openFile('data_new.h5', 'w', 'HiSPARC / KASCADE data')
    hisparc = data.createGroup('/', 'hisparc', 'HiSPARC data')
    kascade = data.createGroup('/', 'kascade', 'KASCADE data')
    coincidences = data.createGroup('/', 'coincidences',
                                    'HiSPARC / KASCADE coincidences')
    data.createTable(hisparc, 'events', HisparcEvent, 'HiSPARC events')
    data.createTable(kascade, 'events', KascadeEvent, 'KASCADE events')
    data.createTable(coincidences, 'events', Coincidence,
                     'Coincidence events')
    data.close()
    print "done."
