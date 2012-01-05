""" Store data from the KASCADE data file using pytables

    This module processes data read from the gzipped data file, provided by
    KASCADE as a courtesy to HiSPARC. The data consists of events taken by
    the KASCADE array, with calculated particle densities at the location
    of our detectors.

    You probably want to use the :func:`read_and_store_data` function.

"""
import gzip
import time

import gpstime

from sapphire.storage import KascadeEvent


class StoreKascadeData():
    def __init__(self, data, hisparc_path, kascade_path, kascade_filename,
                 force=False):

        self.data = data
        self.hisparc = data.getNode(hisparc_path, 'events')

        if kascade_path in data:
            if not force:
                raise RuntimeError("Cancelling simulation; %s already exists?"
                                   % kascade_path)
            else:
                data.removeNode(kascade_path, recursive=True)

        self.kascade = data.createTable(kascade_path, 'events', KascadeEvent,
                                        "KASCADE events", createparents=True)
        self.kascade_filename = kascade_filename

    def read_and_store_data(self):
        """Read and store KASCADE data matching HiSPARC data

        This function looks at the HiSPARC event data in the specified datafile
        and then processes and adds KASCADE data surrounding those events, for
        later coincidence processing.

        """
        # Determine start and end timestamps from HiSPARC data
        try:
            timestamps = self.hisparc.col('timestamp')
            start = gpstime.gps_to_utc(min(timestamps))
            stop = gpstime.gps_to_utc(max(timestamps))
        except IndexError:
            raise RuntimeError("HiSPARC event table is empty")

        print "Processing data from %s to %s" % (time.ctime(start),
                                                 time.ctime(stop))
        self.process_events_in_range(start, stop)

    def process_events_in_range(self, start=None, stop=None):
        """Process KASCADE events in timestamp range

        This function unzips the data file on the fly, reads the data and
        stores it in a pytables table.

        :param start: start of range (timestamp)
        :param stop: end of range (timestamp)

        """
        f = gzip.open(self.kascade_filename)

        while True:
            line = f.readline()
            if not line:
                # no more lines left, EOF
                break

            # break up the line into an array of floats
            data = line.split(' ')
            data = [int(x) for x in data[:4]] + [float(x) for x in data[4:]]

            # KASCADE timestamp
            Gt = data[2]

            # if start and stop are specified, the following boils down to:
            #     start <= Gt < stop
            # but also take start is None and/or stop is None into
            # consideration
            if (start is None or start <= Gt) and (stop is None or Gt < stop):
                self._store_kascade_event(data)
            elif stop is not None and Gt >= stop:
                # timestamp is after explicitly specified stop time, so no need
                # to process the rest of the data
                break

        # flush the table buffers and write them to disk
        self.kascade.flush()

    def _store_kascade_event(self, data):
        tablerow = self.kascade.row

        # read all columns into KASCADE-named variables
        Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
        Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, T200, P200 = data

        tablerow['run_id'] = Irun
        tablerow['event_id'] = Ieve
        tablerow['timestamp'] = Gt
        tablerow['nanoseconds'] = Mmn
        tablerow['ext_timestamp'] = Gt * long(1e9) + Mmn
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

        tablerow.append()
