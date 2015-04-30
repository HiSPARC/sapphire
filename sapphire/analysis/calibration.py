""" Determine calibration values for data

This module can be used to determine calibration values from data.

Determine timing offsets for detectors and stations to correct arrival times.
Determine the PMT response curve to correct the detected number of MIPs.

"""
from ..utils import gauss

from numpy import arange, histogram, percentile, linspace, std, nan, isnan
from scipy.optimize import curve_fit


def determine_detector_timing_offsets(events, station=None):
    """Determine the timing offsets between station detectors.

    :param events: events table of processed events.
    :param station: Station object, to determine number of detectors and
                    relative altitudes.
    :return: list of detector offsets.

    """
    offsets = [nan, nan, nan, nan]
    if not events.nrows:
        return offsets

    t = []
    filters = []
    if station is not None:
        n_detectors = len(station.detectors)
        z = [d.z for d in station.detectors]
    else:
        n_detectors = 4
        z = [0., 0., 0., 0.]

    for id in range(n_detectors):
        t.append(events.col('t%d' % (id + 1)))
        filters.append((events.col('n%d' % (id + 1)) > .05) & (t[id] >= 0.))

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
        offsets[id] = determine_detector_timing_offset(dt, dz)

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
    :return: mean of a gaussian fit to the data corrected for height.

    """
    if not len(dt):
        return nan
    c = .3
    bins = arange(-100 + 1.25, 100, 2.5)
    detector_offset = fit_timing_offset(dt, bins) + dz / c
    return detector_offset


def determine_station_timing_offset(dt, dz=0):
    """Determine the timing offset between station.

    :param dt: a list of time differences between stations (t - t_ref).
    :param dz: height difference between the stations (z - z_ref).
    :return: mean of a gaussian fit to the data corrected for height.

    """
    c = .3
    p = percentile(dt, [2, 98])
    bins = linspace(p[0], p[1], min(int(p[1] - p[0]), 200))
    station_offset = fit_timing_offset(dt, bins) + dz / c
    return station_offset


def fit_timing_offset(dt, bins):
    """Fit the time difference distribution.

    :param dt: a list of time differences between stations (t - t_ref).
    :param bins: bins edges to use for the histogram.
    :return: mean of a gaussian fit to the data corrected for height.

    """
    y, bins = histogram(dt, bins=bins)
    x = (bins[:-1] + bins[1:]) / 2
    try:
        popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0., std(dt)))
        offset = popt[1]
    except RuntimeError:
        offset = nan
    return offset


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
