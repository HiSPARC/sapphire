"""Perform common simulation tasks.

This package contains modules for performing simulations:

:mod:`~sapphire.simulations.base`
    base simulations class that provides the framework for simulations

:mod:`~sapphire.simulations.detector`
    simulate the response of the HiSPARC detectors

:mod:`~sapphire.simulations.groundparticles`
    perform simulations using CORSIKA ground particles as input

:mod:`~sapphire.simulations.ldf`
    perform simulations using lateral distribution functions for particle
    densities

:mod:`~sapphire.simulations.showerfront`
    simple simulations of a shower front

"""
from . import base
from . import detector
from . import groundparticles
from . import ldf
from . import showerfront


__all__ = ['base',
           'detector',
           'groundparticles',
           'ldf',
           'showerfront']
