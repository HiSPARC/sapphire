from struct import unpack

from blocks import (RunHeader, RunEnd, EventHeader, EventEnd,
                    ParticleData, Format, ParticleDataThin, FormatThin)


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
        self.first_particle_index = self._header_index + self.format.subblock_size
        self.last_particle_index = self._end_index - self.format.particle_size

    def get_header(self):
        """Get the Event Header

        :return: an instance of EventHeader

        """
        if not self._header:
            self._header = EventHeader(unpack(self.format.subblock_format,
                                              self._raw_file._contents[self._header_index:
                                                                       self._header_index + self.format.subblock_size]))
        return self._header

    def get_end(self):
        """Get the Event end sub-block

        :return: an instance of EventEnd

        """
        if not self._end:
            self._end = EventEnd(unpack(self.format.subblock_format,
                                                self._raw_file._contents[self._end_index:
                                                                         self._end_index + self.format.subblock_size]))
        return self._end

    def get_particles(self):
        """Generator over particles in the event.

        Use like this::

            for particle in event.get_particles():
                pass

        :yield: each particle in the event

        """
        types = {}
        levels = {}
        done = False
        for sub_block_index in self._raw_file._subblocks_indices(self._header_index,
                                                                 self._end_index):
            for p in range(self.format.particles_per_subblock):
                pos = sub_block_index + p * self.format.particle_size
                particle = self._raw_file._get_particle_record(pos)
                t = int(particle.description / 1000)
                l = particle.description % 10
                if t in types.keys():
                    types[t] += 1
                else:
                    types[t] = 1
                # muon additional information
                if t == 75 or t == 76:
                    continue
                # ignore all observation levels except for nr. 1
                if l != 1:
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
        self._file = open(filename, 'rb')
        self._header_index = None
        self._end_index = None
        self._header = None
        self._end = None
        self._contents = self._file.read()
        self._file.close()
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

        n_blocks = len(self._contents) / self.format.block_size
        for block in range(n_blocks):
            a = unpack('i', self._contents[block * self.format.block_size:block *
                            self.format.block_size + self.format.block_padding_size])[0]
            b = unpack('i', self._contents[(block + 1) * self.format.block_size -
                            self.format.block_padding_size:(block + 1) * self.format.block_size])[0]
            if a != b:
                raise Exception('Block #{block} is not right: ({head}, {tail})'
                                .format(block=block, head=a, tail=b))
        return True

    def get_sub_blocks(self):
        """Get the sub-blocks in the file.

        Normally one would not need this function but it is here because
        I have used it.

        """
        n_blocks = len(self._contents) / self.format.block_size
        for b in xrange(0, n_blocks * self.format.block_size, self.format.block_size):
            for s in xrange(0, self.format.subblocks_per_block):
                pos = b + s * self.format.subblock_size + self.format.block_padding_size
                yield unpack(self.format.subblock_format,
                             self._contents[pos:pos + self.format.subblock_size])

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
            tag = unpack('4s', self._contents[block:block + self.format.field_size])[0]
            if tag == 'EVTH':
                event_head = block
            if tag == 'EVTE':
                yield CorsikaEvent(self, event_head, block)
                event_head = None

    def _subblocks_indices(self, min_sub_block=None, max_sub_block=None):
        """Private method. DO NOT USE! EVER!

        The idea of this method is to get the field indices for the
        beginning and end of the events. It does not unpack the data.

        """
        n_blocks = len(self._contents) / self.format.block_size
        for b in xrange(0, n_blocks * self.format.block_size, self.format.block_size):
            for s in xrange(0, self.format.subblocks_per_block):
                pos = b + s * self.format.subblock_size + self.format.block_padding_size
                if ((not min_sub_block is None and pos <= min_sub_block)
                    or (not max_sub_block is None and pos >= max_sub_block)):
                    continue
                yield pos

    def _get_events(self):
        """Private method. DO NOT USE! EVER!"""
        heads = [b for b in self._subblocks_indices() if unpack('4s', self._contents[b:b + self.format.field_size])[0] == 'EVTH']
        tails = [b for b in self._subblocks_indices() if unpack('4s', self._contents[b:b + self.format.field_size])[0] == 'EVTE']
        return (heads, tails)

    def _get_event_header(self, word):
        """Private method. DO NOT USE! EVER!"""
        return EventHeader(unpack(self.format.subblock_format,
                                  self._contents[word:self.format.subblock_size + word]))

    def _get_event_end(self, word):
        """Private method. DO NOT USE! EVER!"""
        return EventEnd(unpack(self.format.subblock_format,
                                   self._contents[word:self.format.subblock_size + word]))

    def _get_run_indices(self):
        """Get the indices for the start of the run header and end"""
        for block in self._subblocks_indices():
            type = unpack('4s', self._contents[block:block + self.format.field_size])[0]
            if  type == 'RUNH':
                head = block
            elif type == 'RUNE':
                tail = block
        return head, tail

    def _get_run_header(self, word):
        """Get the run header from the contents

        :param word: the index where the run header starts

        """
        return RunHeader(unpack(self.format.subblock_format,
                                self._contents[word:word + self.format.subblock_size]))

    def _get_run_end(self, word):
        """Get the run end from the contents

        :param word: the index where the run end starts

        """
        return RunEnd(unpack(self.format.subblock_format,
                                 self._contents[word:word + self.format.subblock_size]))

    def _get_particle_record(self, word):
        """Private method. DO NOT USE! EVER!"""
        if self.format.particle_format == '7f':
            return ParticleData(unpack(self.format.particle_format,
                                       self._contents[word:word + self.format.particle_size]))
        elif self.format.particle_format == '8f':
            return ParticleDataThin(unpack(self.format.particle_format,
                                           self._contents[word:word + self.format.particle_size]))
        else:
            raise Exception('Unknown particle format: {format}'
                            .format(format=self.format.particle_format))

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
