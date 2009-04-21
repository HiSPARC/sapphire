""" Store data from the KASCADE data file using pytables

    This module processes data read from the gzipped data file, provided by
    KASCADE as a courtesy to HiSPARC. The data consists of events taken by
    the KASCADE array, with calculated particle densities at the location
    of our detectors.

"""
import subprocess

def process_kascade_events(filename, table):
    """Do the actual data processing.

    This function starts a subprocess to unzip the data file, reads the
    data line by line and stores it in a pytables table, row by row.

    Only data starting from July 1st, 14:29 UTC is processed, since that is
    when our detector was turned on.

    Arguments:
    filename    the KASCADE data filename
    table       the destination table

    """
    # start a subprocess which executes gunzip
    p = subprocess.Popen(['gunzip', '-c', filename],
                         stdout=subprocess.PIPE)

    tablerow = table.row

    while True:
        # read a line from the subprocess stdout buffer
        line = p.stdout.readline()
        if not line:
            # no more lines left, EOF
            break

        # break up the line into an array of floats
        data = line.split(' ')
        data = [float(x) for x in data]

        # read all columns into KASCADE-named variables
        Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
        Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, P200, T200 = data

        if Gt >= 1214922585:
            # if timestamp after the HiSPARC experiment started, process
            # the data

            # construct a pytables table row of all the data...
            tablerow['run_id'] = Irun
            tablerow['event_id'] = Ieve
            tablerow['timestamp'] = Gt
            tablerow['nanoseconds'] = Mmn
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

    # flush the table buffers and write them to disk
    table.flush()
