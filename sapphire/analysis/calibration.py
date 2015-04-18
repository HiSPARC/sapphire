""" Determine calibration values for data

This module can be used to determine calibration values from data.

Determine timing offsets for detectors and stations to correct arrival times.
Determine the PMT response curve to correct the detected number of MIPs.

"""
from ..utils import gauss

from numpy import arange, histogram
from scipy.optimize import curve_fit


def determine_detector_timing_offsets(events, station):
    """Determine the timing offsets between station detectors.

    :param events: events table of processed events.
    :param station: Station object.
    :returns: list of detector offsets.

    """
    ref_detector = 1
    t_ref = events.col('t%d' % (ref_detector + 1))
    n_ref = events.col('n%d' % (ref_detector + 1))
    filter = (n_ref > .05) & (t_ref >= 0)
    z = [d.z for d in station.detectors]

    offsets = [0., 0., 0., 0.]
    for detector in range(len(station.detectors)):
        if detector == ref_detector:
            continue
        t = events.col('t%d' % (detector + 1))
        n = events.col('n%d' % (detector + 1))
        dt = (t - t_ref).compress(filter & (n > .05) & (t >= 0))
        dz = z[detector] - z[1]
        offsets[detector] = determine_detector_timing_offset(dt, dz)

    return offsets

def determine_detector_timing_offset(dt, dz=0):
    """Determine the timing offset between station detectors.

    :param dt: a list of time differences between detectors (t - t_ref).
    :param dz: height difference between the detector (z - z_ref).
    :returns: mean of a gaussian fit to the data corrected for height.

    """
    bins = arange(-100 + 1.25, 100, 2.5)
    c = .3
    y, bins = histogram(dt, bins=bins)
    x = (bins[:-1] + bins[1:]) / 2
    try:
        popt, pcov = curve_fit(gauss, x, y, p0=(len(dt), 0., 10.))
        return popt[1] + dz / c
    except (IndexError, RuntimeError):
        return 0.
