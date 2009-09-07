""" Store data from the KASCADE data file using pytables

    This module processes data read from the gzipped data file, provided by
    KASCADE as a courtesy to HiSPARC. The data consists of events taken by
    the KASCADE array, with calculated particle densities at the location
    of our detectors.

"""
import gzip
import time

import gpstime

def process_events(filename, table, start=None, stop=None):
    """Do the actual data processing.

    This function starts a subprocess to unzip the data file, reads the
    data line by line and stores it in a pytables table, row by row.

    Only data starting from July 1st, 14:29 UTC is processed, since that is
    when our detector was turned on.

    Arguments:
    filename    the KASCADE data filename
    table       the destination table

    """
    f = gzip.open(filename)
    tablerow = table.row

    while True:
        # read a line from the subprocess stdout buffer
        line = f.readline()
        if not line:
            # no more lines left, EOF
            break

        # break up the line into an array of floats
        data = line.split(' ')
        data = [float(x) for x in data]

        # read all columns into KASCADE-named variables
        Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
        Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, P200, T200 = data

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
        elif stop is not None and stop < Gt:
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

    Arguments:
    hisparc             HiSPARC event table
    kascade             KASCADE event table
    kascadefile         KASCADE data file

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
        print "data.  This happens, but please check your timestamps:"
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
