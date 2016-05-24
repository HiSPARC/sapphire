""" Determine calibration values for data

This module can be used to determine calibration values from data.

Determine timing offsets for detectors and stations to correct arrival times.
Determine the PMT response curve to correct the detected number of MIPs.

"""
from __future__ import division

from datetime import datetime, timedelta
from itertools import tee, izip, combinations, chain

from numpy import (arange, histogram, percentile, linspace, std, nan, isnan,
                   sqrt, abs, sum)
from scipy.optimize import curve_fit

from ..clusters import HiSPARCStations, HiSPARCNetwork
from ..utils import gauss, round_in_base, memoize, get_active_index, pbar, c
from ..api import Station
from ..transformations.clock import datetime_to_gps, gps_to_datetime


def determine_detector_timing_offsets(events, station=None):
    """Determine the timing offsets between station detectors.

    :param events: events table of processed events.
    :param station: :class:`~sapphire.clusters.Station` object, to determine
        number of detectors and relative altitudes.
    :return: list of detector offsets.

    """
    offsets = [nan, nan, nan, nan]
    if not events.nrows:
        return offsets

    t = []
    filters = []
    if station is not None:
        n_detectors = len(station.detectors)
        station.cluster.set_timestamp(events[0]['timestamp'])
        z = [d.get_coordinates()[2] for d in station.detectors]
    else:
        n_detectors = 4
        z = [0., 0., 0., 0.]

    for id in range(n_detectors):
        t.append(events.col('t%d' % (id + 1)))
        filters.append((events.col('n%d' % (id + 1)) > .3) & (t[id] >= 0.))

    if n_detectors == 2:
        ref_id = 1
    else:
        ref_id = determine_best_reference(filters)

    for id in range(n_detectors):
        if id == ref_id:
            offsets[id] = 0.
            continue
        dt = (t[id] - t[ref_id]).compress(filters[id] & filters[ref_id])
        dz = z[id] - z[ref_id]
        offsets[id], _ = determine_detector_timing_offset(dt, dz)

    # If all except reference are nan, make reference nan.
    if sum(isnan(offsets)) == 3:
        offsets = [nan, nan, nan, nan]

    # Try to make detector 2 the reference point, if it is not nan.
    if not isnan(offsets[1]):
        ref = offsets[1]
        offsets = [o - ref for o in offsets]

    return offsets


def determine_detector_timing_offset(dt, dz=0):
    """Determine the timing offset between station detectors.

    :param dt: a list of time differences between detectors (t - t_ref).
    :param dz: height difference between the detector (z - z_ref).
    :return: mean of a gaussian fit to the data corrected for height, and
             the error of the mean.

    """
    if not len(dt):
        return nan, nan
    p = round_in_base(percentile(dt.compress(abs(dt) < 100), [0.5, 99.5]), 2.5)
    bins = arange(p[0] + 1.25, p[1], 2.5)
    detector_offset, detector_offset_error = fit_timing_offset(dt, bins)
    detector_offset += dz / c
    if abs(detector_offset) > 100:
        detector_offset = nan
    return detector_offset, detector_offset_error


class DetermineStationTimingOffsets(object):
    """Determine the timing offsets between stations"""

    # Maximum distance between station pairs that are included in analysis
    MAX_DISTANCE = 1000  # m
    # Minimum number of timedeltas required to attempt a fit
    MIN_LEN_DT = 100

    def __init__(self, stations=None, data=None, progress=False,
                 force_stale=False):
        """Initialize the class

        :param stations: list of stations for which to determine offsets.
        :param data: the PyTables datafile with timedelta tables.
        :param progress: if true: show progressbar if true.
        :param force_stale: if true: do not get network information from API.

        """
        self.data = data
        self.progress = progress
        self.force_stale = force_stale
        if stations is not None:
            self.cluster = HiSPARCStations(stations, skip_missing=True,
                                           force_stale=self.force_stale)
        else:
            self.cluster = HiSPARCNetwork(force_stale=self.force_stale)

    def read_dt(self, station, ref_station, start, end):
        """Read timedeltas from HDF5 file"""

        table_path = '/time_deltas/station_%d/station_%d' % (ref_station,
                                                             station)
        table = self.data.get_node(table_path, 'time_deltas')
        ts0 = datetime_to_gps(start)  # noqa
        ts1 = datetime_to_gps(end)  # noqa
        return table.read_where('(timestamp >= ts0) & (timestamp < ts1)',
                                field='delta')

    @memoize
    def _get_gps_timestamps(self, station):
        """Get timestamps of station gps changes"""
        return Station(station,
                       force_stale=self.force_stale).gps_locations['timestamp']

    @memoize
    def _get_electronics_timestamps(self, station):
        """Get timestamps of station electronics (hardware) changes"""
        return Station(station,
                       force_stale=self.force_stale).electronics['timestamp']

    def _get_cuts(self, station, ref_station):
        """Get cuts for determination of offsets

        Get a list of events (new gps location, new electronics)
        that (may) cause a large shift in station timing offset
        :param station: station number
        :param ref_station: reference station number
        :return: list of datetime objects

        """
        cuts = {self._datetime(gps_to_datetime(ts))
                for ts in chain(self._get_gps_timestamps(station),
                                self._get_gps_timestamps(ref_station),
                                self._get_electronics_timestamps(station),
                                self._get_electronics_timestamps(ref_station))}
        today = self._datetime(datetime.now())
        cuts = sorted(list(cuts) + [today])
        return cuts

    @memoize
    def _get_r_dz(self, date, station, ref_station):
        """Determine r and dz at date

        :param date: date for which to get the distances.
        :param station,ref_station: station numbers of the station pair.
        :return: tuple containing the horizontal and vertical distances.

        """
        self.cluster.set_timestamp(datetime_to_gps(date))
        r, _, dz = self.cluster.calc_rphiz_for_stations(
            self.cluster.get_station(ref_station).station_id,
            self.cluster.get_station(station).station_id)
        return r, dz

    def _determine_interval(self, r):
        """Determine interval (number of days) in which to fit timedelta's

        :param r: distrance between stations (m).
        :return: number of days in interval.

        """
        # TODO: determine sensible number of days
        return max(int(r ** 1.12 / 10), 7)

    def _get_left_and_right_bounds(self, cuts, date, days):
        """Determine left and right bounds between cuts

        Offsets are determined per day, so intervals are based on days.
        Cuts are excluded. Start date (left side bound) is the day
        after a cut, end date (right side bound) is the day before a cut.
        The last cut (today) is always *included* in the interval,
        as this is not a cut that influences the timing offset.
        Returns datetime objects with hours, min, sec, msec = 0.

        :param cuts: list of datetime objects.
        :param date: datetime (middle of interval).
        :param days: number of days.
        :return: tuple of datetime objects (left bound, right bound).

        """
        left = get_active_index(cuts, self._datetime(date))

        if left == len(cuts) - 1:
            lower_bound = cuts[left - 1]
            upper_bound = cuts[-1]  # include last day (today) in interval
        else:
            right = min(left + 1, len(cuts) - 1)
            lower_bound = cuts[left]
            upper_bound = cuts[right] - timedelta(1)

        step = timedelta(round(days / 2))
        if days >= (upper_bound - lower_bound).days:
            return lower_bound, upper_bound
        elif date + step > upper_bound:
            return upper_bound - 2 * step, upper_bound
        elif date - step < lower_bound:
            return lower_bound, lower_bound + 2 * step
        else:
            return date - step, date + step

    def determine_first_and_last_date(self, date, station, ref_station):
        """
        Determine first and last date to include in determination of
        station offset around date

        :param date: date around which the bounds are to be determined.
        :param station: station number.
        :param ref_station: reference station number.
        :return: start and end date bounds.

        """
        date = self._datetime(date)
        cuts = self._get_cuts(station, ref_station)
        r, dz = self._get_r_dz(date, station, ref_station)
        interval = self._determine_interval(r)

        return self._get_left_and_right_bounds(cuts, date, interval)

    def _datetime(self, date):
        """Ensure date is a datetime object

        :return: a datetime object with h, m, s, ms = 0.

        """
        return datetime(date.year, date.month, date.day)

    def determine_station_timing_offset(self, date, station, ref_station):
        """Determine the timing offset between a station pair at certain date

        :param date: date for which to determine offset as datetime.date.
        :param station: station number.
        :param ref_station: reference station number.
        :return: station offset and reduced chi squared.

        """
        date = self._datetime(date)
        left, right = self.determine_first_and_last_date(date, station,
                                                         ref_station)
        r, dz = self._get_r_dz(date, station, ref_station)
        dt = self.read_dt(station, ref_station, left, right)
        if len(dt) < self.MIN_LEN_DT:
            s_off, rchi2 = nan, nan
        else:
            s_off, rchi2 = determine_station_timing_offset(dt, dz)

        return s_off, rchi2

    def determine_station_timing_offsets(self, station, ref_station,
                                         start=None, end=None):
        """Determine the timing offsets between a station pair

        :param station: station number.
        :param ref_station: reference station number.
        :param start: datetime.date object.
        :param end: datetime.date object.
        :return: list of station offsets as tuple (timestamp, offset, rchi2).

        """
        if start is None:
            cuts = self._get_cuts(station, ref_station)
            start = self._datetime(cuts[0])
        if end is None:
            end = self._datetime(datetime.now())

        offsets = []
        length = (end - start).days
        for date, _ in pbar(datetime_range(start, end), show=self.progress,
                            length=length):
            ts0 = datetime_to_gps(date)
            s_off, rchi2 = self.determine_station_timing_offset(date, station,
                                                                ref_station)
            offsets.append((ts0, s_off, rchi2))
        return offsets

    def determine_station_timing_offsets_for_date(self, date):
        """Determine the timing offsets between a station pair

        :param date: date for which to determine offsets as datetime.date.
        :return: list of station offsets as tuple
                 (station, ref_station, offset, rchi2).

        """
        station_pairs = self.get_station_pairs_within_max_distance(date)
        offsets = []
        for station, ref_station in station_pairs:
            s_off, rchi2 = self.determine_station_timing_offset(date, station,
                                                                ref_station)
            offsets.append((station, ref_station, s_off, rchi2))
        return offsets

    def get_station_pairs_within_max_distance(self, date=None):
        """Iterator that yields stations pairs that are close to each other"""

        if date is not None:
            self.cluster.set_timestamp(datetime_to_gps(date))
        for so1, so2 in combinations(self.cluster.stations, 2):
            s1, s2 = so1.number, so2.number
            r = self.cluster.calc_distance_between_stations(s1, s2)
            if r <= self.MAX_DISTANCE:
                if s1 < s2:
                    yield s1, s2
                else:
                    yield s2, s1


def determine_station_timing_offset(dt, dz=0):
    """Determine the timing offset between stations.

    :param dt: a list of time differences between stations (t - t_ref).
    :param dz: height difference between the stations (z - z_ref).
    :return: mean of a gaussian fit to the data corrected for height, and
             the error of the mean.

    """
    if not len(dt):
        return nan, nan
    p = percentile(dt, [0.5, 99.5])
    # Bins should at least be 1 ns wide, on average at least 4 counts per bin
    # and at most 200 bins.
    bins = linspace(p[0], p[1], min(int(p[1] - p[0]), len(dt) / 4, 200))
    station_offset, station_offset_error = fit_timing_offset(dt, bins)
    station_offset += dz / c
    if abs(station_offset) > 1000:
        return nan, nan
    return station_offset, station_offset_error


def fit_timing_offset(dt, bins):
    """Fit the time difference distribution.

    :param dt: a list of time differences between stations (t - t_ref).
    :param bins: bins edges to use for the histogram.
    :return: mean of a gaussian fit to the data and the error of the mean.

    """
    y, bins = histogram(dt, bins=bins)
    x = (bins[:-1] + bins[1:]) / 2
    sigma = sqrt(y + 1)
    try:
        popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0., std(dt)),
                               sigma=sigma, absolute_sigma=False)
        offset = popt[1]
        width = popt[2]
        offset_error = width / sqrt(sum(y))
    except RuntimeError:
        offset, offset_error = nan, nan
    return offset, offset_error


def determine_best_reference(filters):
    """Find which detector has most events in common with the others

    :param filters: list of filters for each detector, selecting rows
                    where that detector has data.
    :return: index for the detector that has most rows in common with
             the other detectors.

    """
    lengths = []
    ids = range(len(filters))

    for id in ids:
        idx = [j for j in ids if j != id]
        lengths.append(sum(filters[id] & (filters[idx[0]] |
                                          filters[idx[1]] | filters[idx[2]])))
    return lengths.index(max(lengths))


def datetime_range(start, end, step=1):
    """Generator that splits a date range in (almost) equal intervals

    The yielded interval lengths are integer days
    Spreads remaining days over first intervals

    :param start: date instance
    :param end: date instance
    :param step: the integer number of days in each interval
    :return: a tuple of datetime instances for each interval

    """
    interval = (end - start).days

    number_of_steps = interval // step
    if number_of_steps == 0:
        yield start, end
        return

    remainder = interval % step

    chunk_start = start
    for i in range(number_of_steps):
        chunk_end = chunk_start + timedelta(step + min(1, remainder))
        yield chunk_start, chunk_end
        chunk_start = chunk_end
        remainder = max(0, remainder - 1)


def pairwise(iterable):
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""

    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)
