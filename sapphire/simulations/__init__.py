"""Simulations

This package contains modules for performing simulations:

:mod:`~sapphire.simulations.base`
    Base simulations class that provides the framework for simulations

:mod:`~sapphire.simulations.detector`
    Simulate the response of the HiSPARC detectors

:mod:`~sapphire.simulations.groundparticles`
    Perform simulations using CORSIKA ground particles as input

:mod:`~sapphire.simulations.showerfront`
    Simple simulations of a shower front

"""
from . import base
from . import detector
from . import groundparticles
from . import showerfront


__all__ = ['base',
           'detector',
           'groundparticles',
           'showerfront']
