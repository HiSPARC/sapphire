.. include:: subst.inc

Submit CORSIKA jobs to Stoomboot
================================

In order to quickly get a good sample of simulated showers we use the
Nikhef computer cluster Stoomboot to run multiple jobs simultaneously.
For this purpose a modules has been written that will make this easy.
The :mod:`sapphire.corsika.qsub` module can submit as many jobs as you
want with the parameters that you desire.

The function :func:`sapphire.corsika.qsub.multiple_jobs` can do this for
you. The class :class:`sapphire.corsika.qsub.CorsikaBatch` sets up the
environment and starts the jobs. For instance, running 100 showers with
proton primaries of 10e16 eV on the standard Stoomboot queue with the
default |corsika| configuration::

    multiple_jobs(100, 7, 'proton', 'stbcq')

The easiest way to run jobs is to just call the qsub.py script from the shell::

    $ ./qsub.py

Changing the shower settings can be done at the end of this script after::

    if __name__ == '__main__':


.. _CorsikaQsubModule:

Corsika qsub Module
===================

.. automodule:: sapphire.corsika.qsub
   :members:
   :undoc-members:
