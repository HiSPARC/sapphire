Installation
============

In theory, installing a Python package should be as easy as py.  Er, pie.
Using `Pip <http://www.pip-installer.org/>`_.  Like so::

    $ pip install sapphire

Then, the python package called ``sapphire`` would be retrieved from the
internet.  It would have its dependencies listed and ``pip`` would pull
them in and all should be well.  In fact, sapphire *does* have its
dependencies listed and pip *will* pull them in.  It is only then, that
things start go wrong.  Whether you'll experience difficulties depends on
the operating system you're using and previously installed software.
Before I'll go on describing how to install sapphire itself, we will first
install the prerequisites.


Installing the prerequisites
----------------------------

Follow the instructions below for your operating system of choice.


Mac OS X
^^^^^^^^

If you're using Mac OS X, the easiest way to install open source software
(like Python, the program language we're using) including lots and lots of
great packages, is done using `Homebrew
<http://mxcl.github.com/homebrew/>`_.  Please follow the installation
instructions (really easy) and when done, type the following into a
terminal::

    $ brew install python

This will install Python and Pip.

As of this writing, several of the dependencies listed by sapphire do not
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
``apt-get``, but some of them are outdated.  Use ``apt-get`` for python
pacakges at your own discretion (read: risk).  Pip handles python packages
very well, so I'll give some instructions using Pip.

As of this writing, several of the dependencies listed by sapphire do not
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


Installing sapphire
------------------

Like so::

    $ python setup.py install
