"""Process CORSIKA simulation data.

This package contains modules for working with CORSIKA simulations:

:mod:`~sapphire.corsika.blocks`
    classes for the block and subblocks in CORSIKA data

:mod:`~sapphire.corsika.corsika_queries`
    select CORSIKA showers from a corsika overview

:mod:`~sapphire.corsika.generate_corsika_overview`
    generate an overview table of available CORSIKA simulations

:mod:`~sapphire.corsika.particles`
    convert CORSIKA particle codes to common names

:mod:`~sapphire.corsika.qsub_corsika`
    submit CORSIKA simulations to the Nikhef Stoomboot facility

:mod:`~sapphire.corsika.qsub_store_corsika_data`
    submit store CORSIKA jobs to the Nikhef Stoomboot facility

:mod:`~sapphire.corsika.reader`
    read CORSIKA data files into Python

:mod:`~sapphire.corsika.store_corsika_data`
    convert CORSIKA data files to HDF5 files

:mod:`~sapphire.corsika.units`
    convert values in units used by CORSIKA to HiSPARC units

"""
from . import corsika_queries
from . import particles
from . import reader
from . import units


__all__ = ['corsika_queries',
           'particles',
           'reader',
           'units']
