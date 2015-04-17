""" Determine calibration values for data

This module can be used to determine calibration values from data.

Determine timing offsets for detectors and stations to correct arrival times.
Determine the PMT response curve to correct the detected number of MIPs.

"""
from ..utils import gauss

from numpy import arange, histogram
from scipy.optimize import curve_fit


def determine_detector_timing_offset(dt, dz=0):
    """Determine the offset between station detectors.

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
