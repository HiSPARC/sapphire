"""Convert coordinate and time systems.

Easy transformations between different systems.

:mod:`~sapphire.transformations.angles`
    conversion between different angle notations

:mod:`~sapphire.transformations.axes`
    rotation matrices and conversion between coordinate systems

:mod:`~sapphire.transformations.base`
    conversion between decimal and sexagesimal

:mod:`~sapphire.transformations.celestial`
    conversion between celestial coordinate systems, contains both
    legacy and new astropy functions

:mod:`~sapphire.transformations.clock`
    conversion between different time keeping systems

:mod:`~sapphire.transformations.geographic`
    geographic coordinate transformations (e.g. WGS84 to ENU)


"""
from . import angles, axes, base, celestial, clock, geographic

__all__ = ['angles',
           'axes',
           'base',
           'celestial',
           'clock',
           'geographic']
