Example scripts
===============

In this chapter we'll discuss some scripts showing common examples of how to download and analyse data. We'll start with a short discussion of the general layout of such scripts and how to run them interactively. Then, we'll continue with examples for some common tasks.


First, some tips and tricks
---------------------------

The ``not in globals()`` trick
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A simple script downloading data for a single day can look like:

.. literalinclude:: scripts/simple-download.py
   :caption: simple-download.py

This script can be run from the command line like this::

   $ python <name_of_script.py>

Every time you run the script, the old data will be removed and new data will be downloaded. Except that the data is probably unchanged. What a waste of time and bandwidth! If you do *not* run this from the command line you can run into trouble. For example, the second time you run the script from an IPython console, you can run into an exception:

.. code-block:: python
   :emphasize-lines: 9

   >>> %run simple-download.py
   100%|################################################################|Time: 5.22
   >>> %run simple-download.py
   Traceback (most recent call last):
     File "/Users/david/work/HiSPARC/software/sapphire/doc/scripts/simple-download.py", line 10, in <module>
       data = tables.open_file(DATAFILE, 'w')
     File "/Users/david/anaconda/lib/python2.7/site-packages/tables/file.py", line 315, in open_file
       "close it before reopening in write mode." % filename)
   ValueError: The file 'data.h5' is already opened.  Please close it before reopening in write mode.

The second time the script is run, the file is already opened. The script, however, does not know that. You can close the file at the end of your script, but then you don't have access to the file from the interactive console. It is very useful to run the script, keep the file open, and be able to inspect the file and any variables created by the script. The console has access to all variables declared by the script, but by default the script is not allowed to know about any variables defined in the console. In other words, the script *cannot* know about the opened file. This is fixed by running the script with the ``-i`` option::

   >>> %run -i simple-download.py

The script still does not check for opened files. To do this, one can use the ``globals()`` function. This will return a dictionary containing all variables in the global scope. If ``data`` is defined, then ``'data'`` will be present in that dictionary. We can check for that:

.. literalinclude:: scripts/simple-download-with-globals.py
   :caption: simple-download-with-globals.py
   :emphasize-lines: 10-12

Now, the script can be rerun multiple times:

.. code-block:: pycon

   >>> %run -i simple-download-with-globals.py
   100%|################################################################|Time: 5.29
   >>> %run -i simple-download-with-globals.py

The second time the script is run, it does nothing. This *is* useful, because we can include stuff that needs to be rerun *after* the ``if``-statement.


The ``not in data`` trick
^^^^^^^^^^^^^^^^^^^^^^^^^

Every time we run the previous script from the command-line, or start a new IPython console, the script still overwrites the data file with 'new' data. We need something more! In addition to checking if the file is already opened, the script needs to check if the data is already present in the file. If it is, the script should do nothing. It is important to open the file in *append* mode, not in *write* mode. If you do the latter, the file is overwritten the moment you open it. The script becomes:

.. literalinclude:: scripts/simple-download-with-checks.py
   :caption: simple-download-with-checks.py
   :emphasize-lines: 13-14


The ``if __name__ == '__main__'`` trick
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Finally, most scripts of any decent size include a weird-looking line::

   if __name__ == '__main__':
      do_something()

This is very useful, in fact! A small script invariably grows to a larger size. Most probably, you'll write a few functions that are useful elsewhere. If they are, you can import them from a new script. That, however, causes problems if you did not include the lines above. For example, take this script:

.. literalinclude:: scripts/is_useful.py
   :caption: is_useful.py

This script defines a function to calculate the squares of numbers. We use it to print the squared of the numbers 1, 2, and 3. Now, we want to print the square of 16::

   >>> from is_useful import square
   1 2 3
   1 4 9
   >>> print square(16)
   256

What happened? When we import a module, the code *inside* that module is run as well. We don't want that. When we *run* a script, the special variable ``__name__`` is equal to ``'__main__'``. When we *import* a script, the variable is equal to the name of the module, and *not* equal to ``'__main__'``. So, we can check for that:

.. literalinclude:: scripts/is_useful_and_importable.py
   :caption: is_useful_and_importable.py

If we run the script stand-alone, we get the expected result::

   >>> %run is_useful_and_importable.py
   1 2 3
   1 4 9

But we can now also import the script without any unintended side effects::

   >>> from is_useful_and_importable import square
   >>> print square(16)
   256

It is good practice to always include the ``__name__ == '__main__'`` check.


Common tasks
------------

Download event summary data for a few stations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: scripts/download_esd_events.py


Download coincidences for a few stations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: scripts/download_esd_coincidences.py


Download and reconstruct directions for coincidences
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: scripts/download_and_reconstruct_coincidences.py


Plot zenith angle distribution
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. literalinclude:: scripts/plot_zenith_distribution.py
