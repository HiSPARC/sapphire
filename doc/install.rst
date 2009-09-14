.. include:: subst.inc

Installation of the |hisparc| framework
=======================================

The |hisparc| framework is a `Python <http://www.python.org/>`_ package
which only deals with the HiSPARC data and leaves general things like data
storage and plotting to other modules.  Therefore, you'll need to install
Python and several libraries before you can use the framework.
Furthermore, due to the fact that the framework is still in early alpha
and development happens in a quick pace, there are no code releases yet.
Of course, you can grab the sources from the source repository, which is
in Bazaar.

Currently, the instructions are for Ubuntu Linux and Microsoft Windows.
If you want to submit instructions for other platforms, you're very
welcome to do so.  If you find errors or omissions, please inform us so we
can update the documentation!


Prerequisites
-------------

You'll need to install the following in order to make use of the |hisparc|
framework:

`Python 2.6 <http://www.python.org/>`_
    The Python programming language

`NumPy <http://numpy.scipy.org/>`_
    A fundamental package needed for scientific computing with Python

`SciPy <http://scipy.org/>`_
    A software package for mathematics, science, and engineering

`Matplotlib <http://matplotlib.sourceforge.net/>`_
    A Python 2D plotting library which produces publication quality
    figures

`PyTables <http://www.pytables.org/>`_
    A package for managing hierarchical datasets and designed to
    efficiently and easily cope with extremely large amounts of data

`MySQLdb <http://mysql-python.sourceforge.net/>`_
    An interface to the popular MySQL database server, needed to connect
    to our current eventwarehouse

And, preferably:

`IPython <http://ipython.scipy.org/>`_
    An enhanced interactive Python shell, providing a comprehensive
    environment for interactive and exploratory computing

To be able to download and track the framework sources yourself and
optionally take part in development, you'll need:

`Bazaar <http://bazaar-vcs.org/>`_
    a distributed version control system that adapts to the way you want
    to work


Installation in Ubuntu Linux 9.04, the *Jaunty Jackalope*
---------------------------------------------------------

Installation in any Linux distribution is rather painless. Especially in
Debian-based distributions, like Ubuntu, it is as easy as::

    $ sudo apt-get install ipython python-scipy python-matplotlib python-mysqldb python-setuptools
    $ sudo easy_install tables

That's all for the prerequisites. For checking out the framework sources::

    $ bzr checkout sftp://<loginname>@login.nikhef.nl/project/hisparc/bzr/framework/trunk framework


Installation in Microsoft Windows XP
------------------------------------

Installation in Windows has become a *lot* easier now that `Python(x,y)
<http://pythonxy.com/>`_ has become a mature software package.  I'll
gladly accept a nice writeup of the installation procedure.


General usage of the framework in Python scripts
------------------------------------------------

To use the |hisparc| framework from a python script, make sure the source
path is in your python path, like so::

    import sys
    sys.path.append('/path/to/framework')

Or, alternatively, you can add the path to the $PYTHONPATH environment
variable.
