"""Simulation and Analysis Program Package for HiSPARC Research

SAPPHiRE simplifies data access, simulations and analysis for the
`HiSPARC <http://www.hisparc.nl>`_ experiment.  It was born out of a
combination of the old framework and the collection of simulation and
analysis scripts developed by David Fokkema for his PhD thesis work.

The following packages and modules are included:

:mod:`~sapphire.analysis`
    package containing analysis-related modules

:mod:`~sapphire.api`
    publicdb api interface

:mod:`~sapphire.clusters`
    definitions for HiSPARC detectors, stations and clusters

:mod:`~sapphire.corsika`
    package containing CORSIKA simulation related modules

:mod:`~sapphire.gpstime`
    conversion functions for converting GPS time to UTC and vice versa

:mod:`~sapphire.kascade`
    work on KASCADE data

:mod:`~sapphire.publicdb`
    public data access

:mod:`~sapphire.simulations`
    package containing simulation-related modules

:mod:`~sapphire.storage`
    storage-related definitions

:mod:`~sapphire.time_util`
    GPS date/time utility functions

:mod:`~sapphire.transformations`
    geographic coordinate transformations

"""
from . import analysis
from . import api
from . import clusters
from . import corsika
from . import esd
from . import gpstime
from . import kascade
from . import publicdb
from . import simulations
from . import storage
from . import time_util
from . import transformations


__all__ = ['analysis',
           'api',
           'clusters',
           'corsika',
           'esd',
           'gpstime',
           'kascade',
           'publicdb',
           'simulations',
           'storage',
           'time_util',
           'transformations']
