"""
This module provides the functionality to read corsika files. It provides the following main classes:

* :ref:`CorsikaFile`: The file class provides a generator over all events in the file.
* :ref:`CorsikaEvent`: The event class that provides a generator over all particles at ground.

and the following classes that correspond to the sub-blocks defined in
the corsika manual:

* :class:`RunHeader`
* :class:`RunTrailer`
* :class:`EventHeader`
* :class:`EventTrailer`
* :class:`ParticleData`
* :class:`CherenkovData`

Issues
======

This module does not handle platform dependent issues such as byte
ordering (endianness) and field size. This was the result of an
afternoon hack and has only been tested with files generated using 32
bit corsika files on a linux system compiled with gfortran.

* **Field Size**: I'm assuming corsika files were generated in 32 bit mode. In this case, each field in the file is 4 bytes long. If the files were generated in 64 bit mode, each field will be (I think) 64 bits. `NumPy <http://numpy.scipy.org/>`_ has support for 64 bit floats.
* **Endianness**: There is no check for byte ordering. It can be added using Python's `struct module <http://docs.python.org/library/struct.html#struct-alignment>`_.
* **Memory mapping**: This module reads the entire corsika file and stores it in memory. This might be a problem with large files and can be solved by using a memory-mapped file.
* **Thinning**: This module does not handle unthinned showers.

* **Special Particles**: This module currently ignores all special (book-keeping) particles like for muon additional information and history.

More Info
=========
For short information on fortran unformatted binary files, take a look
at http://paulbourke.net/dataformats/fortran/

For detailed information on the corsika format, check the 'Outputs'
chapter in the corsika user manual. In particular, check the 'Normal
Particle Output' section.

Authors
=======

Javier Gonzalez <jgonzalez@ik.fzk.de>
"""

from CorsikaFile import *
import CorsikaBlocks
