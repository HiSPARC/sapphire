from struct import unpack

from CorsikaBlocks import (RunHeader, RunTrailer, EventHeader, EventTrailer,
                           ParticleData, Format, ParticleDataThin, FormatThin)


class CorsikaEvent(object):
    def __init__(self, raw_file, header_index, trailer_index):
        """
        CorsikaEvent constructor.

        The user never calls this. The CorsikaFile does.
        """
        self.fRawFile = raw_file
        self.fHeaderIndex = header_index
        self.fTrailerIndex = trailer_index
        self.format = self.fRawFile.format
        self.fFirstParticle = self.fHeaderIndex + self.format.subblock_size
        self.fLastParticle = self.fTrailerIndex - self.format.particle_size
        self.fHeader = None
        self.fTrailer = None

    def GetHeader(self):
        """
        Returns an instance of EventHeader
        """
        if not self.fHeader:
            self.fHeader = EventHeader(unpack(self.format.subblock_format,
                                              self.fRawFile.fContents[self.fHeaderIndex:self.format.subblock_size + self.fHeaderIndex]))
        return self.fHeader

    def GetTrailer(self):
        """
        Returns an instance of EventTrailer
        """
        if not self.fTrailer:
            self.fTrailer = EventTrailer(unpack(self.format.subblock_format,
                                                self.fRawFile.fContents[self.fTrailerIndex:self.format.subblock_size + self.fTrailerIndex]))
        return self.fTrailer

    def GetParticles(self):
        """
        Generator over particles in the event.

        Use like this::

            for particle in event.GetParticles():
                pass
        """

        types = {}
        levels = {}
        done = False
        for sub_block_index in self.fRawFile._SubBlocksIndices(self.fHeaderIndex,
                                                               self.fTrailerIndex):
            for p in range(self.format.particles_per_subblock):
                pos = sub_block_index + p * self.format.particle_size
                particle = self.fRawFile._GetParticleRecord(pos)
                t = int(particle.fDescription / 1000)
                l = particle.fDescription % 10
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
        """
        String representation (a summary of the event)
        """
        out = self.GetHeader().__str__()
        out += "\n"
        out += self.GetTrailer().__str__()
        return out


class CorsikaFile(object):
    def __init__(self, filename):
        """
        CorsikaFile constructor

        It takes a file name as argument
        """
        self.fFilename = filename
        self.fFile = open(filename, 'rb')
        self.fEvents = None
        self.fRuns = None
        self.fContents = self.fFile.read()
        self.fFile.close()
        self.format = Format()


    def Check(self):
        """
        Some basic sanity checks.

        Fortran unformatted files are written in 'blocks'. Each block
        has a header and trailer. They both contain the same
        information: the number of bytes in the block.

        This function only checks if there is an integer number of
        blocks in the file and if the header and trailer are equal.

        Here would be the place to dynamically check for endiannes and
        field size.
        """
        if len(self.fContents) % self.format.block_size != 0:
            raise Exception('File "{name}" does not have an integer number '
                            'of blocks!'.format(name=self.fFilename))

        n_blocks = len(self.fContents) / self.format.block_size
        for block in range(n_blocks):
            a = unpack('i', self.fContents[block * self.format.block_size:block *
                            self.format.block_size + self.format.block_padding_size])[0]
            b = unpack('i', self.fContents[(block + 1) * self.format.block_size -
                            self.format.block_padding_size:(block + 1) * self.format.block_size])[0]
            if a != b:
                raise Exception('Block #{block} is not right: ({head}, {tail})'
                                .format(block=block, head=a, tail=b))

    def GetSubBlocks(self):
        """
        Get the sub-blocks in the file.

        Normally one would not need this function but it is here because
        I have used it.
        """
        n_blocks = len(self.fContents) / self.format.block_size
        for b in xrange(0, n_blocks * self.format.block_size, self.format.block_size):
            for s in xrange(0, self.format.subblocks_per_block):
                pos = b + s * self.format.subblock_size + self.format.block_padding_size
                yield unpack(self.format.subblock_format,
                             self.fContents[pos:pos + self.format.subblock_size])

    def GetEvents(self):
        """
        Generator of Events

        This method is a generator over the events in the file.
        Use it like this::

            for event in my_file.GetEvents():
                pass
        """
        event_head = None
        for block in self._SubBlocksIndices():
            tag = unpack('4s', self.fContents[block:block + self.format.field_size])[0]
            if tag == 'EVTH':
                event_head = block
            if tag == 'EVTE':
                yield CorsikaEvent(self, event_head, block)
                event_head = None

    def _SubBlocksIndices(self, min_sub_block=None, max_sub_block=None):
        """
        Private method. DO NOT USE! EVER!

        The idea of this method is to get the field indices for the
        beginning and end of the events. It does not unpack the data.
        """
        n_blocks = len(self.fContents) / self.format.block_size
        for b in xrange(0, n_blocks * self.format.block_size, self.format.block_size):
            for s in xrange(0, self.format.subblocks_per_block):
                pos = b + s * self.format.subblock_size + self.format.block_padding_size
                if ((not min_sub_block is None and pos <= min_sub_block)
                    or (not max_sub_block is None and pos >= max_sub_block)):
                    continue
                yield pos

    def _GetEvents(self):
        """
        Private method. DO NOT USE! EVER!
        """
        heads = [b for b in self._SubBlocksIndices() if unpack('4s', self.fContents[b:b + self.format.field_size])[0] == 'EVTH']
        tails = [b for b in self._SubBlocksIndices() if unpack('4s', self.fContents[b:b + self.format.field_size])[0] == 'EVTE']
        return (heads, tails)

    def _GetEventHeader(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return EventHeader(unpack(self.format.subblock_format,
                                  self.fContents[word:self.format.subblock_size + word]))

    def _GetEventTrailer(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return EventTrailer(unpack(self.format.subblock_format,
                                   self.fContents[word:self.format.subblock_size + word]))

    def _GetRunHeader(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return RunHeader(unpack(self.format.subblock_format,
                                self.fContents[word:self.format.subblock_size + word]))

    def _GetRunTrailer(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return RunTrailer(unpack(self.format.subblock_format,
                                 self.fContents[word:self.format.subblock_size + word]))

    def _GetParticleRecord(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        if self.format.particle_format == '7f':
            return ParticleData(unpack(self.format.particle_format,
                                       self.fContents[word:self.format.particle_size + word]))
        elif self.format.particle_format == '8f':
            return ParticleDataThin(unpack(self.format.particle_format,
                                           self.fContents[word:self.format.particle_size + word]))
        else:
            raise Exception('Unknown particle format: {format}'
                            .format(format=self.format.particle_format))

    def Blocks():
        pass


class CorsikaFileThin(CorsikaFile):
    def __init__(self, filename):
        """
        CorsikaFileThin constructor

        It takes a file name as argument
        """
        super(CorsikaFileThin, self).__init__(filename)
        self.format = FormatThin()
