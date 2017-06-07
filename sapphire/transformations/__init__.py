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

:mod:`~sapphire.transformations.celestial_astropy`


"""
from . import angles
from . import axes
from . import base
from . import clock
from . import geographic

"""
    Celestial is an old, deprecated package, it has now been
    replaced with more accurate astropy transformations.
    
    If astropy is present, then celestial_astropy will be used
    and then celestial will be overridden.
    
    If it is not present
    then celestial will be used so that users without astropy
    may still use SAPPHIRE and a warning will be issued.
"""
try:
    import astropy
    from . import celestial_astropy as celestial
except ImportError:
    from . import celestial
    import warnings
    warnings.warn("Astropy not present, falling back on deprecated celestial", ImportWarning)

__all__ = ['angles',
           'axes',
           'base',
           'celestial',
           'clock',
           'geographic']
