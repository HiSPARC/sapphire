"""Data analysis

This package contains modules for performing various analysis tasks:

:mod:`~sapphire.analysis.coincidence_queries`
    Filter coincidences analysed with
    :class:`~sapphire.analysis.coincidences.CoincidenceESD`

:mod:`~sapphire.analysis.coincidences`
    Search for coincidences between HiSPARC stations

:mod:`~sapphire.analysis.core_reconstruction`
    Shower core and shower size reconstruction

:mod:`~sapphire.analysis.direction_reconstruction`
    EAS direction reconstruction

:mod:`~sapphire.analysis.find_mpv`
    Module for finding the MIP peak MPV in pulseintegral data

:mod:`~sapphire.analysis.landau`
    Landau distribution

:mod:`~sapphire.analysis.process_events`
    Process HiSPARC events

:mod:`~sapphire.analysis.reconstructions`
    Perform shower reconstructions

"""
from . import coincidence_queries
from . import coincidences
from . import core_reconstruction
from . import direction_reconstruction
from . import find_mpv
from . import landau
from . import process_events
from . import reconstructions


__all__ = ['coincidence_queries',
           'coincidences',
           'core_reconstruction',
           'direction_reconstruction',
           'find_mpv',
           'landau',
           'process_events',
           'reconstructions']
