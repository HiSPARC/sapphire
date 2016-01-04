SAPPHiRE — A Framework for HiSPARC
===================================

Introduction
------------

.. image:: http://img.shields.io/pypi/v/hisparc-sapphire.svg
   :target: https://pypi.python.org/pypi/hisparc-sapphire/
.. image:: http://img.shields.io/badge/license-GPLv3-blue.svg
   :target: https://github.com/HiSPARC/sapphire/blob/master/LICENSE
.. image:: http://img.shields.io/travis/HiSPARC/sapphire/master.svg
   :target: https://travis-ci.org/HiSPARC/sapphire
.. image:: http://img.shields.io/coveralls/HiSPARC/sapphire/master.svg?label=coveralls
   :target: https://coveralls.io/r/HiSPARC/sapphire
.. image:: http://img.shields.io/codecov/c/github/HiSPARC/sapphire/master.svg?label=codecov
   :target: https://codecov.io/github/HiSPARC/sapphire

SAPPHiRE is a Simulation and Analysis Program Package for `HiSPARC
<http://www.hisparc.nl/>`_ Research and Education.  It was created in the
process of completing the PhD research of David Fokkema.  The history of this
repository contains the complete simulation, analysis and plot generation code
that formed the basis for David's `thesis
<http://www.nikhef.nl/pub/services/biblio/theses_pdf/thesis_D_Fokkema.pdf>`_.
Arne de Laat took over development of SAPPHiRE while working on his own PhD
research.

This repository is created with a sole purpose in mind: to enable HiSPARC
students, teachers and researchers to easily gain access to the data and
perform common simulation and analysis tasks.  Historically, starting work
on the data, or extending an existing analysis code, has involved
elaborate installation instructions, heavy customizations to the software,
countless hours going over opaque parts of code and a general feeling of
anguish and despair.  SAPPHiRE's ultimate goal: no more of that.


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


Version release
---------------

Important: First check if the last commit passes the tests on Travis CI!

To release a new version modify the version number in ``setup.py``. Then
create a commit for the new release with a title like 'Bump version to vX.Y.Z'
and a message that contains a summary of the most important changes since the
last release. Then tag the commit and push it to GitHub::

   $ git tag vX.Y.Z
   $ git push --tags

Then upload the new version to PyPI (this requires the ``wheel`` package)::

   $ python setup.py sdist bdist_wheel upload

The latest version is then available from PyPI.
