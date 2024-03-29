SAPPHiRE — A Framework for HiSPARC
===================================

Introduction
------------

.. image:: https://img.shields.io/pypi/v/hisparc-sapphire
   :target: https://pypi.python.org/pypi/hisparc-sapphire/
.. image:: https://img.shields.io/badge/license-GPLv3-blue
   :target: https://github.com/HiSPARC/sapphire/blob/master/LICENSE
.. image:: https://img.shields.io/github/actions/workflow/status/HiSPARC/sapphire/tests.yml?branch=master
   :target: https://github.com/HiSPARC/sapphire/actions

SAPPHiRE is a Simulation and Analysis Program Package for `HiSPARC
<https://www.hisparc.nl/>`_ Research and Education.  It was created in the
process of completing the PhD research of David Fokkema.  The history of this
repository contains the complete simulation, analysis and plot generation code
that formed the basis for David's `thesis
<https://www.nikhef.nl/pub/services/biblio/theses_pdf/thesis_D_Fokkema.pdf>`_.
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

Required: Python. pip will take care of dependencies, but installing
numpy, scipy and pytables from a python distribution is preferred. We use
miniconda, which includes the conda package manager.

First, `install conda <https://docs.conda.io/en/latest/miniconda.html>`_
and optionally create a virtualenv::

    $ conda create --name hisparc python numpy scipy pytables
    $ source activate hisparc

or alternatively just install the dependencies::

    $ conda install numpy scipy pytables sphinx

Then, using pip::

    $ pip install hisparc-sapphire

This should install sapphire with all requirements. More extensive
installation instructions are available in the documentation in the
``doc/`` directory.  You can compile them using Sphinx, or you can
follow this link: https://docs.hisparc.nl/sapphire/.

To check if it worked start Python and load the package:

.. code-block:: python

    import sapphire

You're done!


Development
-----------

Install python (preferably using conda) as described above but clone
the sapphire repo instead of installing using pip::

    $ git clone https://github.com/HiSPARC/sapphire.git
    $ cd sapphire
    $ pip install -e .[dev]


Version release
---------------

Important: First check if the last commit passes the tests on GitHub Actions!

To release a new version modify the version number in ``setup.py``. Then
create a commit for the new release with a title like 'Bump version to vX.Y.Z'
and a message that contains a summary of the most important changes since the
last release. Then tag the commit and push it to GitHub::

   $ git tag vX.Y.Z
   $ git push --tags

Then upload the new version to PyPI (this requires the ``build``, ``wheel`` and ``twine``
packages)::

   $ python -m build
   $ twine upload dist/hisparc-sapphire-X.Y.Z.tar.gz
   $ twine upload dist/hisparc_sapphire-X.Y.Z-py3-none-any.whl

The latest version is then available from PyPI.
