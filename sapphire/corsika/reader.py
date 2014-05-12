from struct import unpack
import warnings

from blocks import (RunHeader, RunEnd, EventHeader, EventEnd,
                    ParticleData, Format, ParticleDataThin, FormatThin,
                    particle_data)


class CorsikaEvent(object):
    def __init__(self, raw_file, header_index, end_index):
        """CorsikaEvent constructor.

        The user never calls this. The CorsikaFile does.

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
        for sub_block_index in self._raw_file._subblocks_indices(
                self._header_index, self._end_index):
            for p in range(self.format.particles_per_subblock):
                pos = sub_block_index + p * self.format.particle_size
                particle = self._raw_file._get_particle_record_tuple(pos)
                type = particle[6]  # particle type
                level = particle[9]  # observation level

                # skip padding, used to fill a subblock
                if type == 0:
                    continue
                # muon additional information
                if type in [75, 76]:
                    warnings.warn('Ignoring muon additional information.',
                                  UserWarning)
                    continue
                # ignore all observation levels except for nr. 1
                if level != 1:
                    warnings.warn('Only observation level 1 will be read!',
                                  UserWarning)
                    continue

                yield particle

    def __str__(self):
        """String representation (a summary of the event)"""

        out = self.get_header().__str__()
        out += "\n"
        out += self.get_end().__str__()
        return out


class CorsikaFile(object):

    """Corsika output file handler

    This class will probide an interface for Corsika output files.
    Allowing you go get at the events and particles in the file.
    This class is meant for unthinned simulations.

    """

    def __init__(self, filename):
        """CorsikaFile constructor

        :param: the filemame of the CORSIKA data file

        """
        self._filename = filename
        with open(filename, 'rb') as corsika_data:
            self._contents = corsika_data.read()
        self._header_index = None
        self._end_index = None
        self._header = None
        self._end = None
        self.format = Format()

    def check(self):
        """Check DAT file format

        Some basic sanity checks.

        Fortran unformatted files are written in 'blocks'. Each block
        has a header and end. They both contain the same
        information: the number of bytes in the block.

        This function only checks if there is an integer number of
        blocks in the file and if the header and end are equal.

        Here would be the place to dynamically check for endiannes and
        field size.

        """
        if len(self._contents) % self.format.block_size != 0:
            raise Exception('File "{name}" does not have an integer number '
                            'of blocks!'.format(name=self._filename))
        block_size = self.format.block_size
        padding = self.format.block_padding_size
        n_blocks = len(self._contents) / block_size
        for block in range(n_blocks):
            a = unpack('i', self._contents[block * block_size:
                                           block * block_size + padding])[0]
            b = unpack('i', self._contents[(block + 1) * block_size - padding:
                                           (block + 1) * block_size])[0]
            if a != b:
                raise Exception('Block #{block} is not right: ({head}, {tail})'
                                .format(block=block, head=a, tail=b))
        return True

    def get_sub_blocks(self):
        """Get the sub-blocks in the file.

        Normally one would not need this function but it is here because
        I have used it.

        """
        block_size = self.format.block_size
        subblock_size = self.format.subblock_size
        n_blocks = len(self._contents) / block_size
        for b in xrange(0, n_blocks * block_size, block_size):
            for s in xrange(0, self.format.subblocks_per_block):
                pos = b + s * subblock_size + self.format.block_padding_size
                yield unpack(self.format.subblock_format,
                             self._contents[pos:pos + subblock_size])

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
            endblocktag = block + self.format.field_size
            tag = unpack('4s', self._contents[block:endblocktag])[0]
            if tag == 'EVTH':
                event_head = block
            elif tag == 'EVTE':
                yield CorsikaEvent(self, event_head, block)
                event_head = None

    def _subblocks_indices(self, min_sub_block=None, max_sub_block=None):
        """Private method. DO NOT USE! EVER!

        The idea of this method is to get the field indices for the
        beginning and end of the events. It does not unpack the data.

        """
        block_size = self.format.block_size
        subblock_size = self.format.subblock_size
        n_blocks = len(self._contents) / block_size
        for b in xrange(0, n_blocks * block_size, block_size):
            for s in xrange(0, self.format.subblocks_per_block):
                pos = b + s * subblock_size + self.format.block_padding_size
                if ((not min_sub_block is None and pos <= min_sub_block) or
                        (not max_sub_block is None and pos >= max_sub_block)):
                    continue
                yield pos

    def _get_events(self):
        """Get the start and end blocks indices for all events"""

        heads = []
        tails = []
        for block in self._subblocks_indices():
            endblock = block + self.format.field_size
            type = unpack('4s', self._contents[block:endblock])[0]
            if type == 'EVTH':
                heads.append(block)
            elif type == 'EVTE':
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
            endblock = block + self.format.field_size
            type = unpack('4s', self._contents[block:endblock])[0]
            if type == 'RUNH':
                head = block
            elif type == 'RUNE':
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

    def _unpack_subblock(self, word):
        """Unpack a subblock block

        :param word: the index where the subblock starts

        """
        endword = word + self.format.subblock_size
        return unpack(self.format.subblock_format,
                      self._contents[word:endword])

    def _unpack_particle(self, word):
        """Unpack a particle block

        :param word: the index where the particle subblock starts

        """
        endword = word + self.format.particle_size
        return unpack(self.format.particle_format,
                      self._contents[word:endword])

    def Blocks():
        pass


class CorsikaFileThin(CorsikaFile):

    """Corsika thinned output file handler

    Same as the unthinned output handler, but with support for
    the different format, particles also have the weight property.
    This class is meant for thinned simulations.

    """

    def __init__(self, filename):
        """CorsikaFileThin constructor

        It takes a filename as argument

        """
        super(CorsikaFileThin, self).__init__(filename)
        self.format = FormatThin()

    def _get_particle_record(self, word):
        """Get a particle from the contents as a ParticleDataThin instance"""

        return ParticleDataThin(self._unpack_particle(word))
