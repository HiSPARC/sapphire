"""Read CORSIKA data files.

This provides functionality to read CORSIKA output
files with `Python <www.python.org>`_. It provides the following main
classes:

* :class:`~sapphire.corsika.reader.CorsikaFile`: The file class provides a
  generator over all events in the file.
* :class:`~sapphire.corsika.reader.CorsikaEvent`: The event class that
  provides a generator over all particles at ground.

and the following classes that correspond to the sub-blocks defined in
the CORSIKA manual:

* :class:`~sapphire.corsika.blocks.RunHeader`
* :class:`~sapphire.corsika.blocks.RunEnd`
* :class:`~sapphire.corsika.blocks.EventHeader`
* :class:`~sapphire.corsika.blocks.EventEnd`
* :class:`~sapphire.corsika.blocks.ParticleData`
* :class:`~sapphire.corsika.blocks.CherenkovData`

Additionally version for thinned showers are available:

* :class:`CorsikaFileThin`
* :class:`~sapphire.corsika.blocks.ParticleDataThin`
* :class:`~sapphire.corsika.blocks.CherenkovDataThin`


Issues
======

This module does not handle platform dependent issues such as byte
ordering (endianness) and field size. This was the result of an
afternoon hack and has only been tested with files generated using
32 bit CORSIKA files on a linux system compiled with gfortran.

* **Field Size**: According to the CORSIKA user manual section 10.2
  all quantities are written as single precision real numbers
  independently of 32-bit or 64-bit, so each field in the file
  should be 4 bytes long.
* **Endianness**: There is no check for byte ordering. It can be added
  using Python's `struct module
  <http://docs.python.org/library/struct.html#struct-alignment>`_.
* **Special Particles**: This module currently ignores all special
  (book-keeping) particles like for muon additional information and
  history.


More Info
=========

For short information on fortran unformatted binary files, take a look
at http://paulbourke.net/dataformats/reading/

For detailed information on the CORSIKA format, check the 'Outputs'
chapter in the CORSIKA user manual. In particular, check the 'Normal
Particle Output' section.


Authors
=======

- Javier Gonzalez
- Arne de Laat <adelaat@nikhef.nl>

"""

import os
import warnings

from struct import unpack

from .blocks import (
    EventEnd,
    EventHeader,
    Format,
    FormatThin,
    ParticleData,
    ParticleDataThin,
    RunEnd,
    RunHeader,
    particle_data,
    particle_data_thin,
)


class CorsikaEvent:
    def __init__(self, raw_file, header_index, end_index):
        """CorsikaEvent constructor

        The user never calls this. The CorsikaFile does.

        :param raw_file: :class:`CorsikaFile` object.
        :param header_index: index where the event header starts.
        :param end_index: index where the event end starts.

        """
        self._raw_file = raw_file
        self._header_index = header_index
        self._end_index = end_index
        self._header = None
        self._end = None
        self.format = self._raw_file.format
        self.first_particle_index = header_index + self.format.subblock_size
        self.last_particle_index = end_index - self.format.particle_size

    def get_header(self):
        """Get the Event Header

        :return: an instance of EventHeader

        """
        if not self._header:
            self._header = self._raw_file._get_event_header(self._header_index)
        return self._header

    def get_end(self):
        """Get the Event end sub-block

        :return: an instance of EventEnd

        """
        if not self._end:
            self._end = self._raw_file._get_event_end(self._end_index)
        return self._end

    def get_particles(self):
        """Generator over particles in the event.

        .. note::
            This generator filters out additional muon information and
            all particles at observation levels other than 1.

        Use like this::

            for particle in event.get_particles():
                pass

        :yield: each particle in the event

        """
        for sub_block_index in self._raw_file._subblocks_indices(self._header_index, self._end_index):
            for particle in self._raw_file._get_particles(sub_block_index):
                particle_type = particle[6]
                observation_level = particle[9]

                # skip padding, used to fill a subblock
                if particle_type == 0:
                    continue
                # muon additional information
                if particle_type in [75, 76]:
                    warnings.warn('Ignoring muon additional information.')
                    continue
                # ignore all observation levels except for nr. 1
                if observation_level != 1:
                    warnings.warn('Only observation level 1 will be read!')
                    continue

                yield particle

    def __repr__(self):
        return f'{self.__class__.__name__}({self._raw_file!r}, {self._header_index!r}, {self._end_index!r})'


class CorsikaFile:
    """CORSIKA output file handler

    This class will probide an interface for CORSIKA output files.
    Allowing you go get at the events and particles in the file.
    This class is meant for unthinned simulations.

    """

    def __init__(self, filename):
        """CorsikaFile constructor

        :param: the filemame of the CORSIKA data file

        """
        self._filename = filename
        self._size = os.path.getsize(self._filename)
        self._file = open(self._filename, 'rb')
        self._header_index = None
        self._end_index = None
        self._header = None
        self._end = None
        self.format = Format()

    def finish(self):
        """Close the opened CORSIKA data file"""

        self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.finish()

    def check(self):
        """Check DAT file format

        Some basic sanity checks.

        Fortran unformatted files are written in 'blocks'. Each block
        has a header and end. They both contain the same
        information: the number of bytes in the block.

        This function only checks if there is an integer number of
        blocks in the file and if the header and end are equal.

        Here would be the place to dynamically check for endianness and
        field size.

        """
        if self._size % self.format.block_size != 0:
            raise ValueError(f'File "{self._filename}" does not have an integer number of blocks!')
        block_size = self.format.block_size
        padding = self.format.block_padding_size
        n_blocks = self._size // block_size
        for block in range(n_blocks):
            self._file.seek(block * block_size)
            a = unpack('i', self._file.read(padding))[0]
            self._file.seek((block + 1) * block_size - padding)
            b = unpack('i', self._file.read(padding))[0]
            if a != b:
                raise ValueError(f'Block #{block} is not right: ({a}, {b})')
        return True

    def get_sub_blocks(self):
        """Get the sub-blocks in the file.

        Normally one would not need this function but it is here because
        I have used it.

        """
        block_size = self.format.block_size
        subblock_size = self.format.subblock_size
        n_blocks = self._size / block_size
        for b in range(0, n_blocks * block_size, block_size):
            for s in range(self.format.subblocks_per_block):
                pos = b + s * subblock_size + self.format.block_padding_size
                self._file.seek(pos)
                yield unpack(self.format.subblock_format, self._file.read(subblock_size))

    def get_header(self):
        """Get the Run header

        :return: an instance of RunHeader

        """
        if not self._header_index:
            self._header_index, self._end_index = self._get_run_indices()
        if not self._header:
            self._header = self._get_run_header(self._header_index)

        return self._header

    def get_end(self):
        """Get the Run end

        :return: an instance of RunEnd

        """
        if not self._end_index:
            self._header_index, self._end_index = self._get_run_indices()
        if not self._end:
            self._end = self._get_run_end(self._end_index)

        return self._end

    def get_events(self):
        """Generator over the Events in the file

        This method is a generator over the events in the file.
        Use it like this::

            for event in my_file.get_events():
                pass

        """
        event_head = None
        for block in self._subblocks_indices():
            self._file.seek(block)
            tag = unpack('4s', self._file.read(self.format.field_size))[0]
            if tag == b'EVTH':
                event_head = block
            elif tag == b'EVTE':
                yield CorsikaEvent(self, event_head, block)
                event_head = None

    def _subblocks_indices(self, min_sub_block=None, max_sub_block=None):
        """Private method. DO NOT USE! EVER!

        The idea of this method is to get the field indices for the
        beginning and end of the events. It does not unpack the data.

        """
        block_size = self.format.block_size
        subblock_size = self.format.subblock_size
        n_blocks = self._size // block_size
        for b in range(0, n_blocks * block_size, block_size):
            for s in range(self.format.subblocks_per_block):
                pos = b + s * subblock_size + self.format.block_padding_size
                if (min_sub_block is not None and pos <= min_sub_block) or (
                    max_sub_block is not None and pos >= max_sub_block
                ):
                    continue
                yield pos

    def _get_events(self):
        """Get the start and end blocks indices for all events"""

        heads = []
        tails = []
        for block in self._subblocks_indices():
            self._file.seek(block)
            tag = unpack('4s', self._file.read(self.format.field_size))[0]
            if tag == b'EVTH':
                heads.append(block)
            elif tag == b'EVTE':
                tails.append(block)
        return (heads, tails)

    def _get_event_header(self, word):
        """Get the event header from the contents"""

        return EventHeader(self._unpack_subblock(word))

    def _get_event_end(self, word):
        """Get the event header from the contents"""

        return EventEnd(self._unpack_subblock(word))

    def _get_run_indices(self):
        """Get the indices for the start of the run header and end"""
        for block in self._subblocks_indices():
            self._file.seek(block)
            tag = unpack('4s', self._file.read(self.format.field_size))[0]
            if tag == b'RUNH':
                head = block
            elif tag == b'RUNE':
                tail = block
        return head, tail

    def _get_run_header(self, word):
        """Get the run header from the contents"""

        return RunHeader(self._unpack_subblock(word))

    def _get_run_end(self, word):
        """Get the run end from the contents"""

        return RunEnd(self._unpack_subblock(word))

    def _get_particle_record(self, word):
        """Get a particle from the contents as a ParticleData instance"""

        return ParticleData(self._unpack_particle(word))

    def _get_particle_record_tuple(self, word):
        """Get a particle from the contents as a tuple"""

        return particle_data(self._unpack_particle(word))

    def _get_particles(self, word):
        """Get subblock of particles from the contents as tuples"""

        unpacked_particles = self._unpack_particles(word)
        particles = zip(*[iter(unpacked_particles)] * self.format.fields_per_particle)
        return (particle_data(particle) for particle in particles)

    def _unpack_subblock(self, word):
        """Unpack a subblock block

        :param word: the index where the subblock starts

        """
        self._file.seek(word)
        return unpack(self.format.subblock_format, self._file.read(self.format.subblock_size))

    def _unpack_particle(self, word):
        """Unpack a particle block

        :param word: the index where the particle subblock starts

        """
        self._file.seek(word)
        return unpack(self.format.particle_format, self._file.read(self.format.particle_size))

    def _unpack_particles(self, word):
        """Unpack a particles subblock

        :param word: the index where the particle subblock starts

        """
        self._file.seek(word)
        return unpack(self.format.particles_format, self._file.read(self.format.particles_size))

    def __repr__(self):
        return f'{self.__class__.__name__}({self._filename!r})'


class CorsikaFileThin(CorsikaFile):
    """CORSIKA thinned output file handler

    Same as the unthinned output handler, but with support for
    the different format, particles also have the weight property.
    This class is meant for thinned simulations.

    """

    def __init__(self, filename):
        """CorsikaFileThin constructor

        It takes a filename as argument

        """
        super().__init__(filename)
        self.format = FormatThin()

    def _get_particle_record(self, word):
        """Get a particle from the contents as a ParticleDataThin instance"""

        return ParticleDataThin(self._unpack_particle(word))

    def _get_particle_record_tuple(self, word):
        """Get a thinned particle from the contents as a tuple"""

        return particle_data_thin(self._unpack_particle(word))

    def _get_particles(self, word):
        """Get subblock of thinned particles from the contents as tuples"""

        unpacked_particles = self._unpack_particles(word)
        particles = zip(*[iter(unpacked_particles)] * self.format.fields_per_particle)
        return (particle_data_thin(particle) for particle in particles)
