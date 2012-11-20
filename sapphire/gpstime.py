""" GPS time module

    This module provides some conversion functions for converting GPS to
    UTC times and vice versa

    FIXME: This has never been intended as a clean implementation

"""

import time
import calendar


def gps_to_utc(timestamp):
    """Convert GPS time to UTC

    """
    if timestamp < gps_from_string('January 1, 2006'):
        raise Exception, "Dates before January 1, 2006 not implemented!"
    elif timestamp < gps_from_string('January 1, 2009'):
        return timestamp - 14
    elif timestamp < gps_from_string('July 1, 2012'):
        return timestamp - 15
    else:
        return timestamp - 16


def utc_to_gps(timestamp):
    """Convert UTC to GPS time

    """
    if timestamp < utc_from_string('January 1, 2006'):
        raise Exception, "Dates before January 1, 2006 not implemented!"
    elif timestamp < utc_from_string('January 1, 2009'):
        return timestamp + 14
    elif timestamp < utc_from_string('July 1, 2012'):
        return timestamp + 15
    else:
        return timestamp + 16


def utc_from_string(date):
    """Convert a date string to UTC

    """
    t = time.strptime(date, '%B %d, %Y')
    return calendar.timegm(t)


def gps_from_string(date):
    """Convert a date string to GPS time

    """
    t = time.strptime(date, '%B %d, %Y')
    return utc_to_gps(calendar.timegm(t))
