"""
This module provides the functionality to read CORSIKA files.

It provides the following main classes:

* :class:`reader.CorsikaFile`: The file class provides a generator over all events
  in the file.
* :class:`reader.CorsikaEvent`: The event class that provides a generator over all
  particles at ground.

and the following classes that correspond to the sub-blocks defined in
the CORSIKA manual:

* :class:`blocks.RunHeader`
* :class:`blocks.RunTrailer`
* :class:`blocks.EventHeader`
* :class:`blocks.EventTrailer`
* :class:`blocks.ParticleData`
* :class:`blocks.CherenkovData`

Additionally version for thinned showers are available:

* :class:`reader.CorsikaFileThin`
* :class:`blocks.ParticleDataThin`
* :class:`blocks.CherenkovDataThin`


Issues
======

This module does not handle platform dependent issues such as byte
ordering (endianness) and field size. This was the result of an
afternoon hack and has only been tested with files generated using 32
bit CORSIKA files on a linux system compiled with gfortran.

* **Field Size**: According to the Corsika user manual section 10.2
  all quantities are written as single precision real numbers
  independently of 32-bit or 64-bit, so each field in the file
  should be 4 bytes long.
* **Endianness**: There is no check for byte ordering. It can be added using
  Python's `struct module
  <http://docs.python.org/library/struct.html#struct-alignment>`_.
* **Memory mapping**: This module reads the entire CORSIKA file and stores
  it in memory. This might be a problem with large files and can be solved
  by using a memory-mapped file.
* **Thinning**: This module can handle thinned showers with the
  special ..Thin subclasses.
* **Special Particles**: This module currently ignores all special
  (book-keeping) particles like for muon additional information and history.


More Info
=========

For short information on fortran unformatted binary files, take a look
at http://paulbourke.net/dataformats/reading/

For detailed information on the CORSIKA format, check the 'Outputs'
chapter in the CORSIKA user manual. In particular, check the 'Normal
Particle Output' section.

Authors
=======

- Javier Gonzalez <jgonzalez@ik.fzk.de>
- Arne de Laat <adelaat@nikhef.nl>
"""

from reader import *
import blocks
