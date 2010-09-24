""" Store data from the KASCADE data file using pytables

    This module processes data read from the gzipped data file, provided by
    KASCADE as a courtesy to HiSPARC. The data consists of events taken by
    the KASCADE array, with calculated particle densities at the location
    of our detectors.

    You probably want to use the :func:`helper` function, and then move on
    to the :mod:`~hisparc.analysis` package.

"""
import gzip
import time

import gpstime

def process_events(filename, table, start=None, stop=None):
    """Do the actual data processing.

    This function unzips the data file on the fly, reads the data and
    stores it in a pytables table.

    :param filename: the KASCADE data filename
    :param table: the destination table

    """
    f = gzip.open(filename)
    tablerow = table.row

    while True:
        line = f.readline()
        if not line:
            # no more lines left, EOF
            break

        # break up the line into an array of floats
        data = line.split(' ')
        data = [float(x) for x in data]

        # read all columns into KASCADE-named variables
        Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
        Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, T200, P200 = data

        # if start and stop are specified, the following boils down to:
        #     start < Gt < stop
        # but also take start is None and/or stop is None into
        # consideration
        if (start is None or start <= Gt) and (stop is None or Gt < stop):
            # condition matched, so process event

            # construct a pytables table row of all the data...
            tablerow['run_id'] = Irun
            tablerow['event_id'] = Ieve
            tablerow['timestamp'] = Gt
            tablerow['nanoseconds'] = Mmn
            tablerow['ext_timestamp'] = Gt * 1e9 + Mmn
            tablerow['energy'] = EnergyArray
            tablerow['core_pos'] = [Xc, Yc]
            tablerow['zenith'] = Ze
            tablerow['azimuth'] = Az
            tablerow['Num_e'] = Size
            tablerow['Num_mu'] = Nmu
            tablerow['dens_e'] = [He0, He1, He2, He3]
            tablerow['dens_mu'] = [Hmu0, Hmu1, Hmu2, Hmu3]
            tablerow['P200'] = P200
            tablerow['T200'] = T200
            # ...and store it
            tablerow.append()
        elif stop is not None and Gt >= stop:
            # timestamp is after explicitly specified stop time, so no need
            # to process the rest of the data
            break

    # flush the table buffers and write them to disk
    table.flush()

def helper(hisparc, kascade, kascadefile):
    """Helper for quickly processing KASCADE data

    This function looks at the HiSPARC event data in the specified datafile
    and then processes and adds KASCADE data surrounding those events, for
    later coincidence processing.  Also, existing KASCADE data rows are
    inspected to determine the time window.

    :param hisparc: HiSPARC event table
    :param kascade: KASCADE event table
    :param kascadefile: KASCADE data file

    Example::

        >>> import tables
        >>> import hisparc
        >>> data = tables.openFile('kascade.h5', 'a')
        >>> data.createGroup('/', 'kascade', "KASCADE data")
        /kascade (Group) 'KASCADE data'
          children := []
        >>> data.createTable('/kascade', 'events', hisparc.containers.KascadeEvent, "KASCADE events")
        /kascade/events (Table(0,)) 'KASCADE events'
          description := {
          "run_id": Int32Col(shape=(), dflt=0, pos=0),
          "event_id": Int64Col(shape=(), dflt=0, pos=1),
          "timestamp": Time32Col(shape=(), dflt=0, pos=2),
          "nanoseconds": UInt32Col(shape=(), dflt=0, pos=3),
          "ext_timestamp": UInt64Col(shape=(), dflt=0, pos=4),
          "energy": Float64Col(shape=(), dflt=0.0, pos=5),
          "core_pos": Float64Col(shape=(2,), dflt=0.0, pos=6),
          "zenith": Float64Col(shape=(), dflt=0.0, pos=7),
          "azimuth": Float64Col(shape=(), dflt=0.0, pos=8),
          "Num_e": Float64Col(shape=(), dflt=0.0, pos=9),
          "Num_mu": Float64Col(shape=(), dflt=0.0, pos=10),
          "dens_e": Float64Col(shape=(4,), dflt=0.0, pos=11),
          "dens_mu": Float64Col(shape=(4,), dflt=0.0, pos=12),
          "P200": Float64Col(shape=(), dflt=0.0, pos=13),
          "T200": Float64Col(shape=(), dflt=0.0, pos=14)}
          byteorder := 'little'
          chunkshape := (399,)
        >>> hisparc.kascade.helper(data.root.hisparc.cluster_kascade.station_601.events, data.root.kascade.events, 'HiSparc.dat.gz')
        Processing data from Tue Jul  1 16:29:31 2008 to Tue Jul  1 17:15:06 2008

    """
    # Determine start and end timestamps from HiSPARC data
    try:
        hstart = gpstime.gps_to_utc(hisparc[0]['timestamp'])
        hstop = gpstime.gps_to_utc(hisparc[-1]['timestamp'])
    except IndexError:
        print "There is no HiSPARC data yet."
        return
    
    # Determine start and stop timestamps from KASCADE data
    try:
        kstart = kascade[0]['timestamp']
        kstop = kascade[-1]['timestamp']
    except IndexError:
        kstart = None
        kstop = None

    if hstart < kstart and kstart is not None:
        # This should never happen
        print "WARNING: HiSPARC data found which is earlier than KASCADE"
        print "data.  Please check your timestamps:"
        print "First HiSPARC data:", hstart
        print "First KASCADE data:", kstart
        print "The timestamps should not differ by much."

    if kstart is None:
        # There is no KASCADE data yet
        start, stop = hstart, hstop
    else:
        # To avoid duplicates in the case of an aborted previous run,
        # remove rows which share the latest timestamp, because we want to
        # reprocess that timestamp
        l = kascade.getWhereList('timestamp==%d' % kstop)
        kascade.removeRows(start=l[0], stop=l[-1] + 1)
        start, stop = kstop, hstop

    print "Processing data from %s to %s" % (time.ctime(start),
                                             time.ctime(stop))
    process_events(kascadefile, kascade, start, stop)
