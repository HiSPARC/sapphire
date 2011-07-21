import CorsikaBlocks
import struct

from CorsikaBlocks import RunHeader, RunTrailer, EventHeader, EventTrailer, ParticleData, CherenkovData

class CorsikaEvent:
    def __init__(self, raw_file, header_index, trailer_index):
        """
        CorsikaEvent constructor.

        The user never calls this. The CorsikaFile does.
        """
        self.fRawFile = raw_file
        self.fHeaderIndex = header_index
        self.fTrailerIndex = trailer_index
        self.fFirstParticle = header_index + CorsikaBlocks.gSubblockSize
        self.fLastParticle = trailer_index - CorsikaBlocks.gParticleRecordSize
        self.fHeader = None
        self.fTrailer = None

    def GetHeader(self):
        """
        Returns an instance of EventHeader
        """
        if not self.fHeader:
            self.fHeader = CorsikaBlocks.EventHeader(struct.unpack(CorsikaBlocks.gSubblockFormat,
                                                                   self.fRawFile.fContents[self.fHeaderIndex:CorsikaBlocks.gSubblockSize+self.fHeaderIndex]))
        return self.fHeader

    def GetTrailer(self):
        """
        Returns an instance of EventTrailer
        """
        if not self.fTrailer:
            self.fTrailer = CorsikaBlocks.EventTrailer(struct.unpack(CorsikaBlocks.gSubblockFormat,
                                                                     self.fRawFile.fContents[self.fTrailerIndex:CorsikaBlocks.gSubblockSize+self.fTrailerIndex]))
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
        for sub_block_index in self.fRawFile._SubBlocksIndices(self.fHeaderIndex, self.fTrailerIndex):
            for p in range(CorsikaBlocks.gParticlesPerSubblock):
                pos = sub_block_index + p*CorsikaBlocks.gParticleRecordSize
                particle = self.fRawFile._GetParticleRecord(pos)
                t = int(particle.fDescription/1000)
                l = particle.fDescription%10
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


class CorsikaFile:
    def __init__(self, filename):
        """
        CorsikaFile constructor

        It takes a file name as argument
        """
        self.fFilename = filename
        self.fFile =  open(filename, 'rb')
        self.fEvents = None
        self.fRuns = None
        self.fContents = self.fFile.read()
        self.fFile.close()

    #def __del__(self): pass

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
        if len(self.fContents)%CorsikaBlocks.gBlockSize != 0:
            raise Exception('File "{name}" does not have an integer number of blocks!'.format(name=filename))

        n_blocks = len(self.fContents)/CorsikaBlocks.gBlockSize
        for block in range(n_blocks):
            a = struct.unpack('i', self.fContents[block*CorsikaBlocks.gBlockSize:block*CorsikaBlocks.gBlockSize+CorsikaBlocks.gBlockPaddingSize])[0]
            b = struct.unpack('i', self.fContents[(block+1)*CorsikaBlocks.gBlockSize-CorsikaBlocks.gBlockPaddingSize: (block+1)*CorsikaBlocks.gBlockSize])[0]
            if a != b:
                raise Exception('Block #{block} is not right: ({head}, {tail})'.format(block=block, head=a, tail=b))

    def GetSubBlocks(self):
        """
        Get the sub-blocks in the file.

        Normally one would not need this function but it is here because I have used it.
        """
        n_blocks = len(self.fContents)/CorsikaBlocks.gBlockSize
        for b in xrange(0, n_blocks*CorsikaBlocks.gBlockSize, CorsikaBlocks.gBlockSize):
            for s in xrange(0, CorsikaBlocks.gSubBlocksPerBlock):
                pos = b + s*CorsikaBlocks.gSubblockSize + CorsikaBlocks.gBlockPaddingSize
                yield struct.unpack(CorsikaBlocks.gSubblockFormat, self.fContents[pos:pos+CorsikaBlocks.gSubblockSize])

    def GetEvents(self):
        """
        Generator of Events

        This method is a generator over the events in the file. Use it like this::

            for event in my_file.GetEvents():
                pass
        """
        event_head = None
        for block in self._SubBlocksIndices():
            tag = struct.unpack('4s', self.fContents[block:block+CorsikaBlocks.gFieldSize])[0]
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
        n_blocks = len(self.fContents)/CorsikaBlocks.gBlockSize
        for b in xrange(0, n_blocks*CorsikaBlocks.gBlockSize, CorsikaBlocks.gBlockSize):
            for s in xrange(0, CorsikaBlocks.gSubBlocksPerBlock):
                pos = b + s*CorsikaBlocks.gSubblockSize + CorsikaBlocks.gBlockPaddingSize
                if ( not min_sub_block is None and pos <= min_sub_block ) or ( not max_sub_block is None and pos >= max_sub_block ):
                    continue
                yield pos

    def _GetEvents(self):
        """
        Private method. DO NOT USE! EVER!
        """
        heads = [b for b in self._SubBlocksIndices() if struct.unpack('4s', self.fContents[b:b+CorsikaBlocks.gFieldSize])[0] == 'EVTH']
        tails = [b for b in self._SubBlocksIndices() if struct.unpack('4s', self.fContents[b:b+CorsikaBlocks.gFieldSize])[0] == 'EVTE']
        return (heads, tails)


    def _GetEventHeader(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return CorsikaBlocks.EventHeader(struct.unpack(CorsikaBlocks.gSubblockFormat, self.fContents[word:CorsikaBlocks.gSubblockSize+word]))

    def _GetEventTrailer(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return CorsikaBlocks.EventTrailer(struct.unpack(CorsikaBlocks.gSubblockFormat, self.fContents[word:CorsikaBlocks.gSubblockSize+word]))

    def _GetRunHeader(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return CorsikaBlocks.RunHeader(struct.unpack(CorsikaBlocks.gSubblockFormat, self.fContents[word:CorsikaBlocks.gSubblockSize+word]))

    def _GetRunTrailer(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return CorsikaBlocks.RunTrailer(struct.unpack(CorsikaBlocks.gSubblockFormat, self.fContents[word:CorsikaBlocks.gSubblockSize+word]))

    def _GetParticleRecord(self, word):
        """
        Private method. DO NOT USE! EVER!
        """
        return CorsikaBlocks.ParticleData(struct.unpack(CorsikaBlocks.gParticleFormat, self.fContents[word:CorsikaBlocks.gParticleRecordSize+word]))

    def Blocks(): pass
