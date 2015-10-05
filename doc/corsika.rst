.. include:: subst.inc

Corsika Reader Documentation
============================

`CORSIKA <https://web.ikp.kit.edu/corsika/>`_ is an Air Shower
Simulation Program. To be used for detailed simulation of
extensive air showers initiated by high-energy cosmic-ray particles.

This module makes it easy to work with |corsika| output, it was
original written by Javier Gonzalez and adapted for `HiSPARC
<www.hisparc.nl>`_ by Arne de Laat.

The :ref:`CorsikaModule` provides functionality to read |corsika| output
files with `Python <www.python.org>`_. It provides the following main
classes:

* :class:`~sapphire.corsika.reader.CorsikaFile`: The file class provides a
  generator over all events in the file.
* :class:`~sapphire.corsika.reader.CorsikaEvent`: The event class that
  provides a generator over all particles at ground.

This documentation sometimes refers to the |corsika| users manual, this
users manual can be found here `CORSIKA User's Guide
<https://web.ikp.kit.edu/corsika/usersguide/corsika_tech.html>`_

In addition to reading |corsika| output, functionally is provided to
easily submit many |corsika| jobs to the Nikhef batch facility
(Stoomboot).


.. _CorsikaModule:

Corsika Module
--------------

.. automodule:: sapphire.corsika
   :members:


Contents
--------

.. toctree::

   corsika/reader
   corsika/blocks
   corsika/units
   corsika/particles
   corsika/corsika_queries
   corsika/qsub
   corsika/store_corsika_data
