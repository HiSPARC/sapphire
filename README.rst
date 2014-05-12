SAPPHiRE — A Framework for HiSPARC
===================================

Introduction
------------

.. image:: http://img.shields.io/pypi/v/hisparc-sapphire.png
   :target: https://pypi.python.org/pypi/hisparc-sapphire/
.. image:: http://img.shields.io/badge/license-GPLv3-blue.png
   :target: https://github.com/HiSPARC/sapphire/blob/master/LICENSE
.. image:: http://img.shields.io/travis/HiSPARC/sapphire/refactor_simulations.png
   :target: https://travis-ci.org/HiSPARC/sapphire
.. image:: http://img.shields.io/coveralls/HiSPARC/sapphire/refactor_simulations.png
   :target: https://coveralls.io/r/HiSPARC/sapphire?branch=refactor_simulations

SAPPHiRE is a Simulation and Analysis Program Package for `HiSPARC
<http://www.hisparc.nl/>`_ Research and Education.  It was created in the
process of completing the PhD research of David Fokkema.  The history of
this repository contains the complete simulation, analysis and plot
generation code that formed the basis for David's `thesis
<http://www.nikhef.nl/pub/services/biblio/theses_pdf/thesis_D_Fokkema.pdf>`_.

This repository is created with a sole purpose in mind: to enable HiSPARC
students, teachers and researchers to easily gain access to the data and
perform common simulation and analysis tasks.  Historically, starting work
on the data, or extending an existing analysis code, has involved
elaborate installation instructions, heavy customizations to the software,
countless hours going over opaque parts of code and a general feeling of
anguish and despair.  SAPPHiRE's ultimate goal: no more of that.

David has tried very hard to write clean code.  However, as is the nature
of finishing a PhD, severe time constraints prevented him to actually
write well-documented, clean code.  He feels, however, that the presently
available code is a good start.  By releasing it now, it can be used,
accessed, and cleaned up.

In fact, it is probably being cleaned up at this very moment!


Installation
------------

Required: Python with pip, the HDF5 and ATLAS libraries and a
Fortran compiler. 

Then, using pip, simply open a Terminal and do::

    $ pip install hisparc-sapphire

This should install sapphire with all requirements. More extensive
installation instructions are available in the documentation in the
``doc/`` directory.  You can compile them using Sphinx, or you can
follow this link: http://docs.hisparc.nl/sapphire/.

To check if it worked start Python and load the package:

.. code-block:: python

    import sapphire

You're done!
