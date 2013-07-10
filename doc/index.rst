.. include:: subst.inc

Corsika Reader Documentation
============================

`CORSIKA <http://www-ik.fzk.de/~corsika/>`_ is an Air Shower
Simulation Program. To be used for detailed simulation of
extensive air showers initiated by high-energy cosmic-ray particles. 

This module makes it easy to work with |corsika| output, it was
original written by `Javier Gonzalez
<http://www-ik.fzk.de/~jgonzalez/corsika/>`_ and adapted for `HiSPARC
<www.hisparc.nl>`_ by Arne de Laat.

The :ref:`CorsikaModule` provides functionality to read |corsika| output
files with `Python <www.python.org>`_. It provides the following main
classes:

    * :class:`corsika.CorsikaFile.CorsikaFile`: The file class provides a
      generator over all events in the file.
    * :class:`corsika.CorsikaFile.CorsikaEvent`: The event class that provides a
      generator over all particles at ground.

This documentation sometimes refers to the |corsika| users manual, this
users manual can be found here `CORSIKA User's Guide
<http://www-ik.fzk.de/~corsika/usersguide/corsika_tech.html>`_


Installation
------------

To use it, download the `source (zip)
<https://github.com/hisparc/corsika/zipball/master>`_ and uncompress it
somewhere in your PYTHONPATH. Then simply start Python and use
:code:`import corsika` to import the module.
Some examples are included some examples (they require |corsika| output).


Requirements
^^^^^^^^^^^^

- Python 2.7.x
- `NumPy <http://numpy.scipy.org/>`_, which can be installed with
  :code:`pip install numpy`.


Contents
--------

.. toctree::
   :maxdepth: 3
   :numbered:

   corsika


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
