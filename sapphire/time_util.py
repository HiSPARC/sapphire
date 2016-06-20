"""Convert UTC/GPS/local time.

The class in this module is intended to save you from typing trainwrecks
just to convert GPS timestamps to GPS date/time and vice versa.
Trainwrecks invite typos and thus many are easily confused about
UTC/GPS/local time.  No more!

"""
import datetime
import calendar
import time


class GPSTime(object):
    """Date/time utility class."""

    def __init__(self, *args):
        """Instantiate the class.

        The arguments to the __init__ method can take different forms.
        First, it can take year, month, day, hour, minutes, seconds
        arguments.  The year, month and day are mandatory, while hour,
        minutes, seconds default to zero.

        If the argument is a single number, this is interpreted as a
        Unix-like timestamp, in GPS time.

        Examples specifying the exact same time::

            >>> GPSTime(2012, 12, 1)
            >>> GPSTime(1354320000)

        """
        if len(args) == 1:
            self._gpstimestamp = args[0]
        elif len(args) >= 3:
            datetime_ = datetime.datetime(*args)
            timetuple = datetime_.utctimetuple()
            self._gpstimestamp = calendar.timegm(timetuple)
        else:
            raise TypeError("Incorrect arguments")

    def gpstimestamp(self):
        """Return the GPS date/time as a timestamp.

        Example::

            >>> gpstime = GPSTime(2012, 12, 1)
            >>> gpstime.gpstimestamp()
            1354320000

        """
        return self._gpstimestamp

    def description(self):
        """Return the GPS date/time as a string.

        Example::

            >>> gpstime = GPSTime(1354320000)
            >>> gpstime.description()
            'Sat Dec  1 00:00:00 2012'

        """
        timetuple = time.gmtime(self._gpstimestamp)
        return time.asctime(timetuple)

    def datetime(self):
        """Return the GPS date/time as a ``datetime.datetime`` instance.

        Example::

            >>> gpstime = GPSTime(1354320000)
            >>> gpstime.datetime()
            datetime.datetime(2012, 12, 1, 0, 0)

        """
        return datetime.datetime.utcfromtimestamp(self._gpstimestamp)

    def __str__(self):
        """Return sensible description of object."""

        return self.description()
