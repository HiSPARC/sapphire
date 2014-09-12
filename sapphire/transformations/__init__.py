"""Transformation modules

Easy transformations between different systems.

:mod:`~sapphire.angels`
    rotation matrices and conversion betweeen decimal and sexagesimal

:mod:`~sapphire.gpstime`
    conversion functions for converting GPS time to UTC and vice versa

:mod:`~sapphire.geographic`
    geographic coordinate transformations (e.g. WGS84 to ENU)


"""
from . import angles
from . import gpstime
from . import geographic


__all__ = ['angles',
           'gpstime',
           'geographic']
