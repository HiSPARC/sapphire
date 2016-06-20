""" Read and store KASCADE data.

    Read data files provided by the KASCADE collaboration and store them
    in a format compatible with HiSPARC data.

    This module contains the following class:

    :class:`StoreKascadeData`
        Read and store KASCADE data files.

    :class:`KascadeCoincidences`
        Find HiSPARC and KASCADE events that belong together.

"""
import gzip
import time
from os.path import splitext

import numpy as np

from .transformations import clock
from .storage import KascadeEvent


class StoreKascadeData(object):
    def __init__(self, data, kascade_filename, kascade_path='/kascade',
                 hisparc_path=None, force=False, progress=True):
        """Initialize the class.

        :param data: the PyTables datafile
        :param hisparc_path: path to the group containing HiSPARC station data.
        :param kascade_path: path of group where KASCADE data wil be stored.
        :param kascade_filename: filename of the KASCADE data source.
        :param force: overwrite existing KASCADE group if it already exists.
        :param progress: if True, show a progress info will be shown.

        """
        self.data = data
        self.progress = progress

        if hisparc_path is not None:
            self.hisparc = data.get_node(hisparc_path, 'events')
        else:
            self.hisparc = None

        if kascade_path in data:
            if not force:
                raise RuntimeError("Cancelling data storage; %s already exists"
                                   % kascade_path)
            else:
                data.remove_node(kascade_path, recursive=True)

        self.kascade = data.create_table(kascade_path, 'events', KascadeEvent,
                                         "KASCADE events", createparents=True)
        self.kascade_filename = kascade_filename

    def read_and_store_data(self):
        """Read and store KASCADE data matching HiSPARC data

        This function looks at the HiSPARC event data in the specified
        datafile and then processes and adds KASCADE data surrounding
        those events, for later coincidence processing.

        """
        if self.hisparc is not None:
            # Determine start and end timestamps from HiSPARC data
            try:
                timestamps = self.hisparc.col('timestamp')
                start = clock.gps_to_utc(min(timestamps)) - 5
                stop = clock.gps_to_utc(max(timestamps)) + 5
            except IndexError:
                raise RuntimeError("HiSPARC event table is empty")

            if self.progress:
                print "Processing data from %s to %s" % (time.ctime(start),
                                                         time.ctime(stop))
        else:
            start = None
            stop = None
            if self.progress:
                print "Processing all data"

        self._process_events_in_range(start, stop)

    def _process_events_in_range(self, start=None, stop=None):
        """Process KASCADE events in timestamp range

        This function unzips the data file on the fly, reads the data and
        stores it in a pytables table.

        :param start: start of range (timestamp)
        :param stop: end of range (timestamp)

        """
        if splitext(self.kascade_filename)[1] == '.gz':
            f = gzip.open(self.kascade_filename)
        else:
            f = open(self.kascade_filename)

        for line in f:
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
        f.close()

    def _store_kascade_event(self, data):
        """Store a line of KASCADE data in the pytables file

        The stored particle densities are the densities in the plane of the
        shower front. Multiply by cos(zenith) to get the particle density on
        the ground.

        :param data: a list of KASCADE reconstructed shower variables for one
                     event.

        """
        tablerow = self.kascade.row

        # read all columns into KASCADE-named variables
        Irun, Ieve, Gt, Mmn, EnergyArray, Xc, Yc, Ze, Az, Size, Nmu, He0, \
            Hmu0, He1, Hmu1, He2, Hmu2, He3, Hmu3, T200, P200 = data

        tablerow['run_id'] = Irun
        tablerow['event_id'] = Ieve
        tablerow['timestamp'] = Gt
        tablerow['nanoseconds'] = Mmn
        tablerow['ext_timestamp'] = Gt * int(1e9) + Mmn
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


class KascadeCoincidences(object):
    def __init__(self, data, hisparc_group, kascade_group, overwrite=False,
                 ignore_existing=False):
        self.data = data
        self.hisparc_group = data.get_node(hisparc_group)
        self.kascade_group = data.get_node(kascade_group)

        if 'c_index' in self.kascade_group:
            if not overwrite and not ignore_existing:
                raise RuntimeError("I found existing coincidences stored in "
                                   "the KASCADE group")
            elif overwrite:
                data.remove_node(kascade_group, 'c_index')

    def search_coincidences(self, timeshift=0, dtlimit=None, limit=None):
        """Search for coincidences

        This function does the actual searching of coincidences. It uses
        a timeshift to shift the HiSPARC data (we know that these employ
        GPS time, so not taking UTC leap seconds into account). The
        shift will also compensate for delays in the experimental setup.

        The coincidences are stored as an instance variable
        (self.coincidences). Final storage can be done by calling
        :meth:`store_coincidences`.

        :param timeshift: the amount of time the HiSPARC data are shifted (in
            seconds).  Default: 0.
        :param dtlimit: limit on the time difference between hisparc and
            kascade events in seconds.  If this limit is exceeded,
            coincidences are not stored.  Default: None.
        :param limit: limit on the number of KASCADE events investigated.

        """
        h, k = self._get_cached_sorted_id_and_timestamp_arrays()

        # Shift the kascade data instead of the hisparc data. There is less of
        # it, so this is much faster.
        k['ext_timestamp'] += int(-1e9) * timeshift

        if dtlimit:
            # dtlimit in ns
            dtlimit *= 1e9

        coinc_dt, coinc_h_idx, coinc_k_idx = [], [], []

        # First loop through kascade data until we have the first event that
        # occurs _after_ the first hisparc event.
        h_idx = 0
        for k_idx in range(len(k)):
            if k[k_idx][1] > h[h_idx][1]:
                break

        # Limit number of KASCADE events investigated
        if limit:
            max_k_idx = k_idx + limit - 1
        else:
            max_k_idx = np.Inf

        while k_idx <= max_k_idx:
            # Try to get the timestamps of the kascade event and the
            # neighbouring hisparc events.
            try:
                h_t = int(h[h_idx][1])
                k_t = int(k[k_idx][1])
                h_t_next = int(h[h_idx + 1][1])
            except IndexError:
                # Reached beyond the event list.
                break

            # Make sure that while the current hisparc event is _before_ the
            # kascade event, the next hisparc event should occur _after_ the
            # kascade event.  That way, the kascade event is enclosed by
            # hisparc events.
            if k_t > h_t_next:
                h_idx += 1
                continue

            # Calculate the time differences for both neighbors. Make sure to
            # get the sign right. Negative sign: the hisparc event is 'left'.
            # Positive sign: the hisparc event is 'right'.
            dt_left = h_t - k_t
            dt_right = h_t_next - k_t

            # Determine the nearest neighbor and add that to the coincidence
            # list, if dtlimit is not exceeded
            if dtlimit is None or min(abs(dt_left), abs(dt_right)) < dtlimit:
                if abs(dt_left) < abs(dt_right):
                    coinc_dt.append(dt_left)
                    coinc_h_idx.append(h_idx)
                else:
                    coinc_dt.append(dt_right)
                    coinc_h_idx.append(h_idx + 1)
                coinc_k_idx.append(k_idx)

            # Found a match for this kascade event, so continue with the next
            # one.
            k_idx += 1

        self.coincidences = np.rec.fromarrays(
            [coinc_dt, coinc_h_idx, coinc_k_idx], names='dt, h_idx, k_idx')

    def store_coincidences(self):
        self.data.create_table(self.kascade_group, 'c_index',
                               self.coincidences)

    def _get_cached_sorted_id_and_timestamp_arrays(self):
        if not hasattr(self, '_h'):
            self._h = self._get_sorted_id_and_timestamp_array(
                self.hisparc_group)
        if not hasattr(self, '_k'):
            self._k = self._get_sorted_id_and_timestamp_array(
                self.kascade_group)
        return self._h.copy(), self._k.copy()

    def _get_sorted_id_and_timestamp_array(self, group):
        timestamps = group.events.col('ext_timestamp')
        ids = group.events.col('event_id')
        data = np.rec.fromarrays([ids, timestamps],
                                 names='event_id, ext_timestamp')
        data.sort(order='ext_timestamp')
        return data
