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

In a few examples, we will plot the data which we have produced.  We make
use of pylab, which is included with matplotlib.  In OS X and Linux, you
can start an IPython terminal with pylab using::

    $ ipython --pylab

On Windows, where you're probably using Python(x,y), this should be the
default mode.

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


Downloading data
^^^^^^^^^^^^^^^^

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
download data from station 501 and store it in a group with name ``/s501``,
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

As you can see in the reference documentation of
:func:`sapphire.publicdb.download_data`, available by either clicking on
the link or typing in ``help(sapphire.publicdb.download_data)`` in the
interpreter, the function takes six arguments: *file, group, station_id,
start, end* and *get_blobs*.  The last one has the default argument
*False*, and may be omitted, as we have done here.  In our example, we
have opened a file, ``mydata.h5``, and have stored the file *handler* in
the variable ``data``.  So, we passed ``data`` to the function.  The group
name is ``/s501``.  Group names in PyTables are just like folders in a
directory hierarchy.  So, we might have specified
``/path/to/my/hisparc/data/files/for_station/s501``.  It is important to
note that this has absolutely nothing to do with *files*.  Whatever path
you specify, it is all contained *inside* your data file.  Since this is
just a small data file, we have opted for a simple structure.  At the
moment, just one group named ``s501`` at the root of the hierarchy.  Group
names must start with a letter, hence the ``s`` for station.

The *station_id* is simply the station number.  Here, we've chosen to
download data for station 501, located at Nikhef.  The *start* and *end*
parameters specify the date/time range.  Finally, *get_blobs* selects
whether binary data should be downloaded.  This includes the traces of
individual events, error messages, etc.  We've selected the default, which
is *False*.  Binary data is often not neccessary, and will increase the
data size tenfold.


Looking around
^^^^^^^^^^^^^^

If you want to know what groups and tables are contained within the data
file, just ``print`` it::

    >>> print data
    mydata.h5 (File) ''
    Last modif.: 'Sat Dec 29 14:50:55 2012'
    Object Tree:
    / (RootGroup) ''
    /s501 (Group) 'Data group'
    /s501/events (Table(137600,)) 'HiSPARC coincidences table'
    /s501/weather (Table(51513,)) 'HiSPARC weather data'

The *object tree* gives an overview of all groups and tables.  As you can
see, the ``/s501`` group contains two tables, ``events`` and ``weather``.
The events table contains the data from the |hisparc| scintillators, while
the weather table contains data from the (optional) weather station.

To directly access any object in the hierarchy, you can make use of the
``data.root`` object, which points to the root group.  Then, just specify
the remaining path, with dots instead of slashes.  For example, to access
the events table::

    >>> print data.root.s501.events
    /s501/events (Table(137600,)) 'HiSPARC coincidences table'

Of course, we'd like to get some more information.  You can drop the print
statement, and just access the object directly.  PyTables is set up such
that it will give more detailed information whenever you specify the
object directly::

    >>> data.root.s501.events
    /s501/events (Table(137600,)) 'HiSPARC coincidences table'
      description := {
      "event_id": UInt32Col(shape=(), dflt=0, pos=0),
      "timestamp": Time32Col(shape=(), dflt=0, pos=1),
      "nanoseconds": UInt32Col(shape=(), dflt=0, pos=2),
      "ext_timestamp": UInt64Col(shape=(), dflt=0, pos=3),
      "data_reduction": BoolCol(shape=(), dflt=False, pos=4),
      "trigger_pattern": UInt32Col(shape=(), dflt=0, pos=5),
      "baseline": Int16Col(shape=(4,), dflt=-1, pos=6),
      "std_dev": Int16Col(shape=(4,), dflt=-1, pos=7),
      "n_peaks": Int16Col(shape=(4,), dflt=-1, pos=8),
      "pulseheights": Int16Col(shape=(4,), dflt=-1, pos=9),
      "integrals": Int32Col(shape=(4,), dflt=-1, pos=10),
      "traces": Int32Col(shape=(4,), dflt=-1, pos=11),
      "event_rate": Float32Col(shape=(), dflt=0.0, pos=12)}
      byteorder := 'little'
      chunkshape := (704,)

There you go!  But what does it all *mean*?  Well, if you want to get to
the bottom of it, read the `PyTables documentation
<http://pytables.github.com/usersguide/index.html>`_.  We'll give a quick
overview here.

First, this table contains 137600 rows.  In total, there are thirteen
columns: ``event_id``, ``timestamp``, ``nanoseconds``, ``ext_timestamp``,
``data_reduction``, ``trigger_pattern``, ``baseline``, ``std_dev``,
``n_peaks``, ``pulseheights``, ``integrals``, ``traces`` and
``event_rate``.

Each event has a unique [#event_id]_ identifier, ``event_id``.  Each event
has a Unix timestamp in GPS time, *not* UTC.  A `Unix timestamp
<http://en.wikipedia.org/wiki/Unix_time>`_ is the number of seconds that
have passed since January 1, 1970.  The sub-second part of the timestamp
is given in ``nanoseconds``.  The ``ext_timestamp`` is the full timestamp
in ns.  Since there cannot exist another event with the same timestamp,
this field in combination with the station number uniquely identifies the
event.  The ``data_reduction`` flag signifies whether the full PMT trace
(*no* reduction) has been stored, or just the PMT pulse (*reduced*, or
*zero suppression*).  The ``trigger_pattern`` is a binary value containing
the exact trigger condition at the time of the event.  The ``baseline``,
``std_dev``, ``n_peaks``, ``pulseheights`` and ``integrals`` fields are
values derived from the PMT traces.  Each field contains four values, one
for each detector.  If a station only has two detectors, the last two
values for each field are -1.  If the baseline cannot be determined, all
these values are -999.  The ``event_rate`` is the trigger rate at the time
of the event.

We'll get to work with this data in a moment.  First, we'll take a look at
the weather table::

    >>> data.root.s501.weather
    /s501/weather (Table(51513,)) 'HiSPARC weather data'
      description := {
      "event_id": UInt32Col(shape=(), dflt=0, pos=0),
      "timestamp": Time32Col(shape=(), dflt=0, pos=1),
      "temp_inside": Float32Col(shape=(), dflt=0.0, pos=2),
      "temp_outside": Float32Col(shape=(), dflt=0.0, pos=3),
      "humidity_inside": Int16Col(shape=(), dflt=0, pos=4),
      "humidity_outside": Int16Col(shape=(), dflt=0, pos=5),
      "barometer": Float32Col(shape=(), dflt=0.0, pos=6),
      "wind_dir": Int16Col(shape=(), dflt=0, pos=7),
      "wind_speed": Int16Col(shape=(), dflt=0, pos=8),
      "solar_rad": Int16Col(shape=(), dflt=0, pos=9),
      "uv": Int16Col(shape=(), dflt=0, pos=10),
      "evapotranspiration": Float32Col(shape=(), dflt=0.0, pos=11),
      "rain_rate": Float32Col(shape=(), dflt=0.0, pos=12),
      "heat_index": Int16Col(shape=(), dflt=0, pos=13),
      "dew_point": Float32Col(shape=(), dflt=0.0, pos=14),
      "wind_chill": Float32Col(shape=(), dflt=0.0, pos=15)}
      byteorder := 'little'
      chunkshape := (1310,)

We'll let these column names speak for themselves.


Accessing the data
^^^^^^^^^^^^^^^^^^

We can access the data in several ways.  We can address the complete
table, or just one or several rows from it.  We can read out a single
column, or select data based on a query.  Before we do any of that, we'll
save us some typing::

    >>> events = data.root.s501.events

Now, we have stored a short-hand reference to the events table. Let's get
the first event::

    >>> events[0]
    (1L, 1354320004, 670634817L, 1354320004670634817L, True, 196608L,
    [195, 197, 197, 197], [658, 763, 646, 771], [1, 0, 1, 0], [173, 3,
    407, 3], [1603, 0, 4019, 0], [0, 1, 2, 3], 0.9111111164093018)

That's the first event!  It is not, however, immediately clear what
numbers correspond to which columns.  They are in order, however, so you
could find out.  It is often easier to specify the column you're
interested in::

    >>> events[0]['pulseheights']
    array([173,   3, 407,   3], dtype=int16)
    >>> events[0]['n_peaks']
    array([1, 0, 1, 0], dtype=int16)

Which gives us the pulseheights of the first event.  The pulseheights are
16-bit integers (that's the ``dtype=int16``) and are determined after
digitizing the events using an analog-digital converter (ADC).  Each unit
corresponds to about -0.57 mV.  You can tell that the first and third
detectors had relatively large pulseheights, and that they were registered
as a significant *peak* in the signal (see the ``n_peaks`` column).

If you're interested in the pulseheights of *all* events, the fastest way
to do it is to make use of the :meth:`Table.col`  method of the table::

    >>> events.col('pulseheights')
    array([[173,   3, 407,   3],
           [  1,   2, 313, 756],
           [211,   2,   2, 268],
           ...,
           [  2,   2, 529, 327],
           [  7,   2, 318, 249],
           [  3, 250,  14, 416]], dtype=int16)

It is also possible to select the data based on a query.  For example, to
select all events between timestamps 1354320000 and 1354323600 (a one-hour
time span)::

    >>> t0 = 1354320000
    >>> t1 = t0 + 3600
    >>> sel_events = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')
    >>> len(sel_events)
    2836

Thus, we have selected 2836 events.  The variable ``sel_events`` no longer
points to a table.  We can no longer make use of the :meth:`Table.col`
method, but we *can* access all pulseheights in the following way::

    >>> sel_events['pulseheights']
    array([[173,   3, 407,   3],
           [  1,   2, 313, 756],
           [211,   2,   2, 268],
           ...,
           [328, 556, 255,  15],
           [  3, 325, 309,   3],
           [  2,   2, 271, 381]], dtype=int16)

This notation is possible for arrays, but not for tables.  So, for tables,
use the :meth:`Table.col` method.  For arrays, use this special notation.


The :mod:`sapphire.time_util` module
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

|sapphire| includes a handy module: :mod:`sapphire.time_util`.  This saves
you from the hassle of converting timestamps and confusing local time and
GPS (or UTC) time.  It is important to realize that the |hisparc| station
uses a GPS clock, and thus saves all timestamps in GPS time, which is
certainly not local time!  You can look up `GPS time
<http://en.wikipedia.org/wiki/Global_Positioning_System#Timekeeping>`_,
but suffice it to say that it is *almost* equal to UTC time.  The
difference is the leap seconds introduced after 1980.  In January 2013,
GPS time is ahead of UTC by 16 seconds.  We will not reference UTC or
local time, but instead always reference GPS time!

While you tell |sapphire| to download data using year, month, day, hour,
minutes, seconds notation, the events table contains timestamps.  It is
often hard to convert between the two.  For example, the correct
conversion between a GPS date/time and a GPS timestamp is given by::

    >>> calendar.timegm(datetime.datetime(2012, 12, 1, 0).utctimetuple())
    1354320000

Which is, at best, cumbersome.  It is easy, however, to screw up and
inadvertently convert to local time.  For your benefit, we have included
the :class:`sapphire.time_util.GPSTime` class.  You instantiate the class
by giving it a GPS time to work with.  It can either be in date/time
notation, or as a timestamp.  For example, the exact same result is
obtained by these two lines of code::

    >>> sapphire.time_util.GPSTime(2012, 12, 1)
    <sapphire.time_util.GPSTime instance at 0x10b8021b8>
    >>> sapphire.time_util.GPSTime(1354320000)
    <sapphire.time_util.GPSTime instance at 0x10b802128>

If you store the instance, you can then call several methods to convert
the date/time to whatever you require::

    >>> gpstime = sapphire.time_util.GPSTime(2012, 12, 1)
    >>> gpstime.datetime()
    datetime.datetime(2012, 12, 1, 0, 0)
    >>> gpstime.description()
    'Sat Dec  1 00:00:00 2012'
    >>> gpstime.gpstimestamp()
    1354320000

Or, indeed::

    >>> gpstime = sapphire.time_util.GPSTime(1354320000)
    >>> gpstime.datetime()
    datetime.datetime(2012, 12, 1, 0, 0)
    >>> gpstime.description()
    'Sat Dec  1 00:00:00 2012'
    >>> gpstime.gpstimestamp()
    1354320000

It is now easy to select events occuring between 12:00 and 13:00 hours GPS
time on December 2, 2012::

    >>> t0 = sapphire.time_util.GPSTime(2012, 12, 2, 12).gpstimestamp()
    >>> t1 = sapphire.time_util.GPSTime(2012, 12, 2, 13).gpstimestamp()
    >>> t0, t1
    (1354449600, 1354453200)
    >>> sel_events = events.readWhere('(t0 <= timestamp) & (timestamp < t1)')
    >>> len(sel_events)
    2817

And when did the tenth event occur?  Well::

    >>> t = sel_events[9]['timestamp']
    >>> t
    1354449620
    >>> sapphire.time_util.GPSTime(t).description()
    'Sun Dec  2 12:00:20 2012'

To shorten the typing somewhat, remember that in Python you can import
classes from modules into your own namespace::

    >>> from sapphire.time_util import GPSTime
    >>> GPSTime(t).description()
    'Sun Dec  2 12:00:20 2012'


Plotting data
^^^^^^^^^^^^^

Now that we can access the data, we want to visualize it.  Plotting data
is a great way to do that.  And of course, the venerable histogram is
still very useful to condense thousands of events into one display.  Pylab
contains an easy function to do just that: :func:`hist`.  Let's try to
recreate a few graphs as seen on the |hisparc| `data display
<http://data.hisparc.nl/django/show/stations/501/2012/12/1/>`_::

    >>> ph = events.col('pulseheights')
    >>> hist(ph)
    ([array([114794,  17565,   3062,   1009,    502,    285,    170,     86,
               72,     55]), array([112813,  19028,   3339,   1246,    540,    295,    163,    100,
               66,     10]), array([109162,  21833,   4246,   1345,    579,    290,    113,     32,
                0,      0]), array([109996,  21283,   4028,   1285,    581,    251,    133,     43,
                0,      0])], array([  1.00000000e+00,   3.91100000e+02,   7.81200000e+02,
             1.17130000e+03,   1.56140000e+03,   1.95150000e+03,
             2.34160000e+03,   2.73170000e+03,   3.12180000e+03,
             3.51190000e+03,   3.90200000e+03]), <a list of 4 Lists of Patches objects>)

.. image:: images/tutorial-hist-simple.png
   :width: 500px

This will not do.  Firstly, data from the four detectors is pictured as
four side-by-side colored bars.  Secondly, the number of bins is very low;
it is only ten.  Thirdly, the data range continues up to very high values
with hardly any events.

To fix this, we'll make use of several arguments that can be passed to the
:func:`hist` function.  We'll also make use of some NumPy (`documentation
<http://docs.scipy.org/doc/numpy/user/>`_) functionality.  For Pylab
documentation, see the `Matplotlib site <http://matplotlib.org>`_
(`function reference
<http://matplotlib.org/1.2.0/api/pyplot_summary.html>`_).  Try this::

    >>> bins = arange(0, 2001, 20)
    >>> hist(ph, bins, histtype='step', log=True)
    >>> xlabel("Pulseheight [ADC]")
    >>> ylabel("Counts")
    >>> title("Pulseheight histogram (log scale)")

The :func:`arange` function returns an array with range from 0 to 2001 in
steps of 20.  It is necessary to say 2001 (or 2002, 2000.1 or whatever)
and not 2000, if you want the range to be inclusive.  The :func:`hist`
function will then plot a *stepped* histogram with a log scale.  Finally,
we add some labels and a title.  This is the result:

.. image:: images/tutorial-hist-better.png
   :width: 500px



.. rubric:: Footnotes

.. [#event_id]

    Unique in this table.  When data is downloaded for analysis and
    combined with other data into one table, the ``event_id`` will be
    different.  To uniquely define an event, use a station number /
    ``ext_timestamp`` combination.
