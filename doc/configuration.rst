.. include:: subst.inc

Configuration
*************

Various modules (:mod:`sapphire.api`, :mod:`sapphire.esd`, and :mod:`sapphire.publicdb`)
connect to the |hisparc| Public Database to download data, metadata, configurations, etc.
Normally you will want to use the official |hisparc| Public Database, which is hosted
at `http://data.hisparc.nl`. However, during development, tests or other occasions you
may want to use a different server, for example when running a local version of the
publicdb. In that case you can overwrite the base url used. This overwrite is done
using an environment variable: ``PUBLICDB_BASE``. This can be easily set when starting
your Python environment, for example::

    $ PUBLICDB_BASE=http://localhost:8000 python

Or while already running Python::

    from os import environ
    environ['PUBLICDB_BASE'] = 'http://localhost:8000'
