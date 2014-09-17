"""Transformation modules

Easy transformations between different systems.

:mod:`~sapphire.axes`
    rotation matrices and conversion between coordinate systems

:mod:`~sapphire.base`
    conversion betweeen decimal and sexagesimal

:mod:`~sapphire.gpstime`
    conversion functions for converting GPS time to UTC and vice versa

:mod:`~sapphire.geographic`
    geographic coordinate transformations (e.g. WGS84 to ENU)


"""
from . import axes
from . import base
from . import gpstime
from . import geographic


__all__ = ['axes',
           'base',
           'gpstime',
           'geographic']
