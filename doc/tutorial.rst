.. include:: subst.inc

Tutorial
========

SAPPHiRE simplifies data access, simulations and analysis for the `HiSPARC
<http://www.hisparc.nl>`_ experiment.  In this tutorial, we'll try to give
you a feeling for the things you can *do* with |sapphire|.  How can you
download data?  How can you analyze this data?  How can you produce
pulseheight histograms?  How can you calculate the direction of cosmic
rays?  This tutorial will only give you an overview of what's possible
with |sapphire|.  For details on all available classes and methods, please
see the :doc:`sapphire`.

.. note::
    We'll require you to know some basic Python.  If not, you can look at
    the official `Python tutorial <http://docs.python.org/2/tutorial/>`_.
    If you know how to start a Python interpreter, you can probably start
    best with chapter 3 and ignore chapter 2.  Chapters 3 through 6 should
    give you some working knowledge of Python.


First steps
-----------

Whenever you see something like::

    >>> print 'Hello, world'
    Hello, world

it means we're typing in Python code.  The ``>>>`` is the Python *prompt*.
It tells you that the Python interpreter is waiting for you to give it
some instructions.  In the above example, we've typed in ``print 'Hello,
world'``.  The interpreter subsequently executed the code and printed
*Hello, world* on the screen.  Some further examples (try them, if you
like!)::

    >>> a = 2 + 2
    >>> 3 * a
    12

The first thing we'll have to do to start using |sapphire| is to *import*
the |sapphire| module::

    >>> import sapphire

There will be no output if everything is succesful.  It is easy to get
some help from inside the Python terminal.  Just say::

    >>> help(sapphire)

and you'll be presented with a basic help screen.  Instead of
``sapphire``, you can throw in modules, packages, functions, classes,
objects, etc.  Everything in Python has *some* help text associated with
it.  Not all of it is very helpful to a newcomer, hence this tutorial.
All help text is also available in the :doc:`sapphire`.
