"""Simulation and Analysis Program Package for HiSPARC Research

SAPPHiRE simplifies data access, simulations and analysis for the `HiSPARC
<http://www.hisparc.nl>`_ experiment.  It was born out of a combination of the
old framework and the collection of simulation and analysis scripts developed
by David Fokkema for his PhD thesis work.  Development is ongoing, while Arne
de Laat is working on his PhD research.

The following packages and modules are included:

:mod:`~sapphire.analysis`
    package containing analysis-related modules

:mod:`~sapphire.api`
    publicdb api interface

:mod:`~sapphire.clusters`
    definitions for HiSPARC detectors, stations and clusters

:mod:`~sapphire.corsika`
    package containing CORSIKA simulation related modules

:mod:`~sapphire.data`
    scripts for updating local data

:mod:`~sapphire.esd`
    event summary data access

:mod:`~sapphire.kascade`
    work on KASCADE data

:mod:`~sapphire.publicdb`
    public data access

:mod:`~sapphire.qsub`
    submit jobs to Stoomboot

:mod:`~sapphire.simulations`
    package containing simulation-related modules

:mod:`~sapphire.storage`
    storage-related definitions

:mod:`~sapphire.tests`
    code tests

:mod:`~sapphire.time_util`
    GPS date/time utility functions

:mod:`~sapphire.transformations`
    transformations between different systems

:mod:`~sapphire.utils`
    commonly used functions such as a progressbar

"""
from . import (
    analysis,
    api,
    clusters,
    corsika,
    data,
    esd,
    kascade,
    publicdb,
    qsub,
    simulations,
    storage,
    time_util,
    transformations,
    utils,
)
from .analysis.calibration import DetermineStationTimingOffsets, determine_detector_timing_offsets
from .analysis.coincidence_queries import CoincidenceQuery
from .analysis.coincidences import Coincidences, CoincidencesESD
from .analysis.find_mpv import FindMostProbableValueInSpectrum
from .analysis.process_events import (
    ProcessEvents,
    ProcessEventsFromSource,
    ProcessEventsFromSourceWithTriggerOffset,
    ProcessSingles,
    ProcessSinglesFromSource,
    ProcessWeather,
    ProcessWeatherFromSource,
)
from .analysis.process_traces import DataReduction, MeanFilter, TraceObservables
from .analysis.reconstructions import (
    ReconstructESDCoincidences,
    ReconstructESDEvents,
    ReconstructESDEventsFromSource,
)
from .analysis.time_deltas import ProcessTimeDeltas
from .api import Network, Station
from .clusters import HiSPARCNetwork, HiSPARCStations, ScienceParkCluster
from .corsika.corsika_queries import CorsikaQuery
from .esd import download_coincidences, download_data, download_lightning, load_data, quick_download
from .simulations.groundparticles import GroundParticlesSimulation, MultipleGroundParticlesSimulation
from .simulations.ldf import KascadeLdfSimulation, NkgLdfSimulation
from .simulations.showerfront import ConeFrontSimulation, FlatFrontSimulation
from .tests import run_tests
from .transformations.celestial import zenithazimuth_to_equatorial
from .transformations.clock import datetime_to_gps, gps_to_datetime
from .version import __version__  # noqa

__all__ = ['analysis',
           'api',
           'clusters',
           'corsika',
           'data',
           'esd',
           'kascade',
           'publicdb',
           'qsub',
           'simulations',
           'storage',
           'time_util',
           'transformations',
           'utils',
           'determine_detector_timing_offsets',
           'DetermineStationTimingOffsets',
           'CoincidenceQuery',
           'Coincidences', 'CoincidencesESD',
           'FindMostProbableValueInSpectrum',
           'ProcessEvents', 'ProcessEventsFromSource',
           'ProcessEventsFromSourceWithTriggerOffset',
           'ProcessWeather', 'ProcessWeatherFromSource',
           'ProcessSingles', 'ProcessSinglesFromSource',
           'TraceObservables', 'MeanFilter', 'DataReduction',
           'ReconstructESDEvents', 'ReconstructESDEventsFromSource',
           'ReconstructESDCoincidences',
           'ProcessTimeDeltas',
           'Network', 'Station',
           'HiSPARCStations', 'HiSPARCNetwork', 'ScienceParkCluster',
           'CorsikaQuery',
           'quick_download', 'load_data', 'download_data',
           'download_lightning', 'download_coincidences',
           'GroundParticlesSimulation', 'MultipleGroundParticlesSimulation',
           'KascadeLdfSimulation', 'NkgLdfSimulation',
           'FlatFrontSimulation', 'ConeFrontSimulation',
           'run_tests',
           'zenithazimuth_to_equatorial',
           'gps_to_datetime', 'datetime_to_gps'
           ]
