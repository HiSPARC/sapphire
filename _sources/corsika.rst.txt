.. include:: subst.inc

CORSIKA simulations
===================

`CORSIKA <https://web.ikp.kit.edu/corsika/>`_ is an Air Shower
Simulation Program. To be used for detailed simulation of
extensive air showers initiated by high-energy cosmic-ray particles.

This module makes it easy to work with CORSIKA output, it was
original written by Javier Gonzalez and adapted for `HiSPARC
<www.hisparc.nl>`_ by Arne de Laat.

This documentation sometimes refers to the CORSIKA users manual, this
users manual can be found here `CORSIKA User's Guide
<https://www.ikp.kit.edu/corsika/70.php>`_

In addition to reading CORSIKA output, functionally is provided to
easily submit many CORSIKA jobs to the Nikhef batch facility
(Stoomboot).


CORSIKA Module
--------------

.. automodule:: sapphire.corsika
   :members:


Contents
--------

.. toctree::
   :hidden:

   corsika/blocks
   corsika/corsika_queries
   corsika/generate_corsika_overview
   corsika/particles
   corsika/qsub_corsika
   corsika/qsub_store_corsika_data
   corsika/reader
   corsika/store_corsika_data
   corsika/units
