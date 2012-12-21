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


Downloading and accessing |hisparc| data
----------------------------------------

The |sapphire| package comprises multiple modules and packages.  To access
data from the public database, we'll have to import the ``publicdb``
module.  Also, we need the ``tables`` module to actually store the data.
`PyTables <pytables.org>`_ is based on the open HDF5 data format, which is
used by `many (research) institutes
<http://www.hdfgroup.org/HDF5/users5.html>`_.  For example, it is used by
the KNMI and by NASA.  To specify the date and time for which to download
the data, we need the ``datetime`` module.  Thus, we have::

    >>> import tables
    >>> import datetime
    >>> import sapphire.publicdb

Creating an empty data file, with the name ``mydata.h5``, is done easily::

    >>> data = tables.openFile('mydata.h5', 'w')

The ``'w'`` means *write*, which creates a file for writing (and reading).
Mind that this will create an empty file.  If there already was a file
with that name, it will be overwritten!  Alternatively, you can say
``'a'``, which means *append*, thus adding to an existing file without
overwriting its contents.  Finally, you can specify ``'r'`` for
*read-only*.

To download data, we have to specify the date/time *range*.  If we want to
download data from the December 1, 2012 all through December 2, 2012,
we can specify this by typing::

    >>> start = datetime.datetime(2012, 12, 1)
    >>> end = datetime.datetime(2012, 12, 3)

Mind that if we do not specify the hour of day, it is taken to be 00:00
hours.  Thus, there is no data included from December 3.  Alternatively,
we can download data from a two hour interval on December 14 by specifying
the hour of day::

    >>> start = datetime.datetime(2012, 12, 14, 19)
    >>> end = datetime.datetime(2012, 12, 14, 21)

which is from 19:00 to 21:00 hours.  It is important to realize that the
we use a GPS clock, which equal to UTC (up to leap seconds).  So, if we
download data for a station in the Netherlands, we have just said from
20:00 to 22:00 local time.  You can specify the time up to the seconds.

We have not actually done anything yet.  We have just stored our time
window in two arbitrarily-named variables, ``start`` and ``end``.  To
download data from station 501 and store it in a group with name ``s501``,
we can use the :func:`sapphire.publicdb.download_data` function::

    >>> sapphire.publicdb.download_data(data, '/s501', 501, start, end)
    INFO:hisparc.publicdb:2012-12-01 00:00:00 None
    INFO:hisparc.publicdb:Getting server data URL (2012-12-01 00:00:00)
    INFO:hisparc.publicdb:Downloading data...
    INFO:hisparc.publicdb:Storing data...
    INFO:hisparc.publicdb:Done.
    INFO:hisparc.publicdb:2012-12-02 00:00:00 None
    INFO:hisparc.publicdb:Getting server data URL (2012-12-02 00:00:00)
    INFO:hisparc.publicdb:Downloading data...
    INFO:hisparc.publicdb:Storing data...
    INFO:hisparc.publicdb:Done.
