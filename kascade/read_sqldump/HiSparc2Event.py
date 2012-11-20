""" Process HiSPARC messages from a buffer
    This module processes messages from buffer database and
    gets out all available data. This data is stored in a data which
    can then be uploaded to the eventwarehouse.
"""

import struct
import base64
from Event import Event
import EventExportValues


class HiSparc2Event(object, Event):
    def __init__(self, message):
        """ Initialization
            First, determine message type from the argument. Then, check if
            this might be a legacy message. Proceed to unpack the message.
        """

        # invoke constructor of parent class
        Event.__init__(self)

        # get the message field in the message table
        self.message = message[1]

    #--------------------------End of __init__--------------------------#

    def unpackMessage(self):
        pass

    #--------------------------End of unpackMessage--------------------------#

    def parseMessage(self):
        self.unpackMessage()

        # get all event data necessary for an upload.
        self.export_values = EventExportValues.export_values[self.uploadCode]

        return self.getEventData()

    #--------------------------End of parseMessage--------------------------#

    def getEventData(self):
        """    Get all event data necessary for an upload.
            This function parses the export_values variable declared in the EventExportValues
            and figures out what data to collect for an
            upload to the eventwarehouse. It returns a list of
            dictionaries, one for each data element.
        """

        eventdata = []
        for value in self.export_values:
            is_calculated = value[0]
            data_uploadcode = value[1]

            try:
                data = self.__getattribute__(value[2])
            except AttributeError:
                #if not self.version == 21:
                    # This is not a legacy message. Therefore, it should
                    # contain all exported variables, but alas, it
                    # apparently doesn't.
                    #print 'I missed this variable: ', value[2]
                continue

            if data_uploadcode in ['TR1', 'TR2', 'TR3', 'TR4']:
                # encode compressed binary traces for transport over http.
                # The pickled data stream is ascii-safe, but the binary
                # compressed traces are not. At the server side, the
                # blobvalues are base64-decoded.
                data = base64.b64encode(data)

            eventdata.append({
                "calculated": is_calculated,
                "data_uploadcode": data_uploadcode,
                "data": data,
            })

        return eventdata

    #--------------------------End of getEventData--------------------------#

    def unpackSeqMessage(self, fmt=None):
        """    Sequentially unpack message with a format
            This method is used to read from the same buffer multiple times,
            sequentially. A private variable will keep track of the current
            offset. This is more convenient than keeping track of it yourself
            multiple times, or hardcoding offsets.
        """
        if not fmt:
            # This is an initialization call
            self._struct_offset = 0
            return

        if fmt == 'LVstring':
            # Request for a labview string. That is, first a long for the
            # length, then the string itself.
            length, = self.unpackSeqMessage('>L')
            fmt = ">%ds" % length

        # For debugging, keeping track of trailing bytes
        #print len(self.message[self._struct_offset:]), struct.calcsize(fmt)

        data = struct.unpack_from(fmt, self.message,
                                  offset=self._struct_offset)
        self._struct_offset += struct.calcsize(fmt)

        return data

    #--------------------------End of unpackSeqMessage--------------------------#
