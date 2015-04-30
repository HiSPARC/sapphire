.. include:: ../subst.inc

Submit CORSIKA jobs to Stoomboot
================================

In order to quickly get a good sample of simulated showers we use the
Nikhef computer cluster Stoomboot to run multiple jobs simultaneously.
For this purpose a script has been written that will make this easy.
The :mod:`~sapphire.corsika.qsub_corsika` script can submit as many
jobs as you want with the parameters that you desire.

The syntax for calling the script can be seen by calling its help::

    qsub_corsika --help

For example, running 100 showers with proton primaries of 10e16 eV
coming in at 22.5 degrees on the standard Stoomboot queue with the
default |corsika| configuration::

    qsub_corsika 100 16 proton 22.5


Corsika on Stoomboot
====================

.. automodule:: sapphire.corsika.qsub_corsika
   :members:
   :undoc-members:
