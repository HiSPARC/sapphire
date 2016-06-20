"""Perform common data analysis tasks on HiSPARC data.

This package contains modules for performing various analysis tasks:

:mod:`~sapphire.analysis.calibration`
    determine calibration values

:mod:`~sapphire.analysis.coincidence_queries`
    filter coincidences analysed with
    :class:`~sapphire.analysis.coincidences.CoincidencesESD`

:mod:`~sapphire.analysis.coincidences`
    search for coincidences between HiSPARC stations

:mod:`~sapphire.analysis.core_reconstruction`
    shower core and shower size reconstruction

:mod:`~sapphire.analysis.direction_reconstruction`
    EAS direction reconstruction

:mod:`~sapphire.analysis.event_utils`
    get data from processed events

:mod:`~sapphire.analysis.find_mpv`
    find the MIP peak MPV in pulseintegral data

:mod:`~sapphire.analysis.landau`
    Landau distribution

:mod:`~sapphire.analysis.process_events`
    process HiSPARC events

:mod:`~sapphire.analysis.process_traces`
    process HiSPARC traces

:mod:`~sapphire.analysis.reconstructions`
    perform shower reconstructions

:mod:`~sapphire.analysis.time_deltas`
    determine time deltas for station pairs

"""
from . import calibration
from . import coincidence_queries
from . import coincidences
from . import core_reconstruction
from . import direction_reconstruction
from . import event_utils
from . import find_mpv
from . import landau
from . import process_events
from . import process_traces
from . import reconstructions
from . import time_deltas


__all__ = ['calibration',
           'coincidence_queries',
           'coincidences',
           'core_reconstruction',
           'direction_reconstruction',
           'event_utils',
           'find_mpv',
           'landau',
           'process_events',
           'process_traces',
           'reconstructions',
           'time_deltas']
