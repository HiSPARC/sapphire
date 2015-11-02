"""Convert coordinate and time systems.

Easy transformations between different systems.

:mod:`~sapphire.transformations.angles`
    conversion between different angle notations

:mod:`~sapphire.transformations.axes`
    rotation matrices and conversion between coordinate systems

:mod:`~sapphire.transformations.base`
    conversion betweeen decimal and sexagesimal

:mod:`~sapphire.transformations.celestial`
    conversion betweeen celestial coordinate systems

:mod:`~sapphire.transformations.clock`
    conversion between different time keeping systems

:mod:`~sapphire.transformations.geographic`
    geographic coordinate transformations (e.g. WGS84 to ENU)

"""
from . import angles
from . import axes
from . import base
from . import celestial
from . import clock
from . import geographic


__all__ = ['angles',
           'axes',
           'base',
           'celestial',
           'clock',
           'geographic']
