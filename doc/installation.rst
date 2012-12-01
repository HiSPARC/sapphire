.. include:: subst.inc

Installation
============

In theory, installing a Python package should be as easy as py.  Er, pie.
Using `Pip <http://www.pip-installer.org/>`_.  Like so::

    $ pip install sapphire

Then, the python package called ``sapphire`` would be retrieved from the
internet.  It would have its dependencies listed and ``pip`` would pull
them in and all should be well.  In fact, |sapphire| *does* have its
dependencies listed and pip *will* pull them in.  It is only then, that
things start go wrong.  Whether you'll experience difficulties depends on
the operating system you're using and any previously installed software.
But don't worry, we've got you covered.  Before I'll go on describing how
to install |sapphire| itself, we will first install the prerequisites.


Installing the prerequisites
----------------------------

Follow the instructions below for your operating system of choice.


Mac OS X
^^^^^^^^

If you're using Mac OS X, the easiest way to install open source software
(like Python, the program language we're using) including lots and lots of
great packages, is done using `Homebrew
<https://mxcl.github.com/homebrew/>`_.  Please follow the installation
instructions (really easy) and when done, type the following into a
terminal::

    $ brew install python

This will install Python and Pip.

As of this writing, several of the dependencies listed by |sapphire| do not
have their own dependencies listed in a way that pip (or other tools, for
that matter) know how to handle.  Furthermore, matplotlib needs to be
installed all by itself, and its dependencies must be installed before it.

After trial and error, this is the magic incantation which works on my
system::

    $ pip install numpy
    $ pip install numexpr
    $ pip install cython
    $ pip install matplotlib

There are more prerequisites to be installed, but they are correctly
handled by the python package management software.


Windows
^^^^^^^

While Python follows a *batteries included* philosophy with an extensive
standard library, the unix philosophy favors minimalism.  Therefore, it is
custom to install Python and then continue to install additional packages
as needed.  Not so with Windows.  Windows does not have a package manager
which can handle many software packages.  Instead, it favors installing
relatively few software packages encompassing lots of functionality.
Installing a package requires running a separate installer.

Luckily, the `Python(x,y) <http://code.google.com/p/pythonxy/>`_ project
solves this problem by installing not only Python, but a complete
scientific environment.  This means that all are prerequisites are
automatically covered.  Please install Python(x,y).


Debian and derivatives (like Ubuntu)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Debian has a first-class package manager which other distributions and
operating systems have a hard time competing with.  In my personal
opinion, anyway.  Lots of python packages can be installed using
``apt-get``, but some of them might be outdated, depending on the age of
your distribution.  Use ``apt-get`` for python pacakges at your own
discretion (read: risk).  Pip handles python packages very well, so I'll
give some instructions using Pip.

As of this writing, several of the dependencies listed by |sapphire| do not
have their own dependencies listed in a way that pip (or other tools, for
that matter) know how to handle.  Furthermore, matplotlib needs to be
installed all by itself, and its dependencies must be installed before it.

You need to install some development libraries::

    $ sudo apt-get install gfortran
    $ sudo apt-get install python-dev
    $ sudo apt-get install libfreetype6-dev
    $ sudo apt-get install libpng-dev
    $ sudo apt-get install libatlas-base-dev
    $ sudo apt-get install libhdf5-dev

After trial and error, this is the magic incantation which works on my
system::

    $ sudo pip install numpy
    $ sudo pip install numexpr
    $ sudo pip install cython
    $ sudo pip install matplotlib

There are more prerequisites to be installed, but they are correctly
handled by the python package management software.


Installing |sapphire|
---------------------

When |sapphire| reaches 1.0, we will upload it to `PyPI
<http://pypi.python.org>`_ so that Pip knows where to find it.  That is
currently, however, not the case.

There are two scenarios for installing |sapphire|: with or without
fetching the code.


Just let me get to work!
^^^^^^^^^^^^^^^^^^^^^^^^

This scenario does not involve fetching the code.  It will just install
|sapphire|, so that you can get to work.  Quickly, issue:

    $ pip install https://github.com/hisparc/sapphire/zipball/master

Done.  Now get to work.


Let me see the code!
^^^^^^^^^^^^^^^^^^^^

You want to see the code so that you can change it, or follow the progress
of |sapphire|.  If you're interested in the development of |sapphire|, you
can either go to the `GitHub page <https://github.com/hisparc/sapphire/>`_
or install the version control system (we use `Git <http://git-scm.com>`_
yourself.  For that, please see the `GitHub Help pages
<https://help.github.com/articles/set-up-git>`_.

To just download the code and install |sapphire|, first go to
https://github.com/hisparc/sapphire/.  Then, click on the *Zip* button
(see image below).  This will start a download of all the code bundled in
a zip file.

.. image:: images/github-zipball.png

Uncompress the zip file, open a terminal and navigate to the top-level
directory containing the code.  Then issue::

    $ python setup.py install

This takes care of installing the missing dependencies and |sapphire|
itself.


Checking that |sapphire| is installed correctly
-----------------------------------------------

First off, the following is not an exhaustive check.  But it will tell you
if |sapphire| is, in fact, installed on your system and that Python knows
how to find it.

Start a Python session.  You can use a launcher of some type (e.g. the one
from Python(x,y)), or open a terminal and type::

    $ python

Or, if you prefer, `IPython <http://ipython.org>`_::

    $ ipython

In fact, we recommend using IPython for interactive use.  Then, try to
import |sapphire|::

    >>> import sapphire

If this returns without an error message, all is well and |sapphire| is
correctly installed.

.. note:: When you run this check from *inside the |sapphire| code
          directory*, it will always return successfully.  The reason for
          this is that Python also checks the current working directory
          for packges.  So, run this check from e.g. your home directory.
