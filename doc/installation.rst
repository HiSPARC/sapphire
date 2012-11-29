Installation
============

In theory, installing a Python package should be as easy as py.  Er, pie.
Using Pip.  Like so::

    $ pip install sapphire

Then, the python package called ``sapphire`` would be retrieved from the
internet.  It would have its dependencies listed and ``pip`` would pull
them in and all should be well.  In fact, sapphire *does* have its
dependencies listed and pip *will* pull them in.  It is only then, that
things start go wrong.  Whether you'll experience difficulties depends on
the operating system you're using and previously installed software.
Before I'll go on describing how to install sapphire itself, we will first
install the prerequisites.


Mac OS X
--------

Homebrew + pip


Windows
-------

Python(x,y)


Debian / Ubuntu
---------------

apt-get + pip


Using pip
---------

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


Installing sapphire
------------------

Like so::

    $ python setup.py install
