import tables

class HisparcEvent(tables.IsDescription):
    event_id = tables.UInt64Col(pos=0)
    timestamp = tables.Time32Col(pos=1)
    nanoseconds = tables.UInt32Col(pos=2)
    ext_timestamp = tables.UInt64Col(pos=3)
    pulseheights = tables.Int16Col(shape=4, dflt=-9999, pos=4)
    integrals = tables.Int32Col(shape=4, dflt=-9999, pos=5)

class KascadeEvent(tables.IsDescription):
    run_id = tables.IntCol(pos=0)
    event_id = tables.Int64Col(pos=1)
    timestamp = tables.Time32Col(pos=2)
    nanoseconds = tables.UInt32Col(pos=3)
    ext_timestamp = tables.UInt64Col(pos=4)
    energy = tables.FloatCol(pos=5)
    core_pos = tables.FloatCol(shape=2, pos=6)
    zenith = tables.FloatCol(pos=7)
    azimuth = tables.FloatCol(pos=8)
    Num_e = tables.FloatCol(pos=9)
    Num_mu = tables.FloatCol(pos=10)
    dens_e = tables.FloatCol(shape=4, pos=11)
    dens_mu = tables.FloatCol(shape=4, pos=12)
    P200 = tables.FloatCol(pos=13)
    T200 = tables.FloatCol(pos=14)

class Coincidence(tables.IsDescription):
    hisparc_event_id = tables.UInt64Col(pos=0)
    kascade_event_id = tables.UInt64Col(pos=1)
    hisparc_timestamp = tables.Time32Col(pos=2)
    hisparc_nanoseconds = tables.UInt32Col(pos=3)
    hisparc_ext_timestamp = tables.UInt64Col(pos=4)
    hisparc_pulseheights = tables.Int16Col(shape=4, dflt=-9999, pos=5)
    hisparc_integrals = tables.Int32Col(shape=4, dflt=-9999, pos=6)
    kascade_timestamp = tables.Time32Col(pos=7)
    kascade_nanoseconds = tables.UInt32Col(pos=8)
    kascade_ext_timestamp = tables.UInt64Col(pos=9)
    kascade_energy = tables.FloatCol(pos=10)
    kascade_core_pos = tables.FloatCol(shape=2, pos=11)
    kascade_zenith = tables.FloatCol(pos=12)
    kascade_azimuth = tables.FloatCol(pos=13)
    kascade_Num_e = tables.FloatCol(pos=14)
    kascade_Num_mu = tables.FloatCol(pos=15)
    kascade_dens_e = tables.FloatCol(shape=4, pos=16)
    kascade_dens_mu = tables.FloatCol(shape=4, pos=17)
    kascade_P200 = tables.FloatCol(pos=18)
    kascade_T200 = tables.FloatCol(pos=19)

def create_tables():
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


if __name__ == '__main__':
    create_tables()
