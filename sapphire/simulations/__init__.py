"""Simulations

This package contains modules for performing simulations:

:mod:`~sapphire.simulations.base`
    Base simulations class that provides the framework for simulations

:mod:`~sapphire.simulations.groundparticles`
    Perform simulations using CORSIKA ground particles as input

"""
from . import base
from . import groundparticles


__all__ = ['base',
           'groundparticles']
