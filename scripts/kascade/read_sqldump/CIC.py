"""
    Process HiSPARC messages from a buffer.
    This module processes the CIC event message
"""

__author__ = "thevinh"
__date__ = "$17-sep-2009"

import struct
import datetime
from zlib import compress

from HiSparc2Event import HiSparc2Event
from legacy import unpack_legacy_message
import EventExportValues


class CIC(HiSparc2Event):
    def __init__(self, message):
        """ Initialization
            Proceed to unpack the message.
        """
        # invoke constructor of parent class
        HiSparc2Event.__init__(self, message)

        # init the trigger rate attribute
        self.eventrate = 0

                self.uploadCode = 'CIC'

    #--------------------------End of __init__--------------------------#

    def parseMessage(self):
        # get database flags
        tmp = struct.unpack("B", self.message[0:1])[0]
        if tmp <= 3:
            unpack_legacy_message(self)
        else:
            self.unpackMessage()

        # get all event data necessary for an upload.
        self.export_values = EventExportValues.export_values[self.uploadCode]

        return self.getEventData()

    #--------------------------End of parseMessage--------------------------#

    def unpackMessage(self):
        """    Unpack a buffer message
            This routine unpacks a buffer message written by the LabVIEW DAQ
            software version 3.0 and above. Version 2.1.1 doesn't use a version
            identifier in the message. By including one, we can account for
            different message formats.

            Hopefully, this code is cleaner and thus easier to understand than
            the legacy code. However, you'll always have to be careful with the
            format strings.
        """

        # Initialize sequential reading mode
        self.unpackSeqMessage()

        self.version, self.database_id, self.data_reduction, \
        self.eventrate, self.num_devices, self.length, \
        gps_second, gps_minute, gps_hour, gps_day, gps_month, gps_year, \
        self.nanoseconds, self.time_delta, self.trigger_pattern = \
        self.unpackSeqMessage('>2BBfBH5BH3L')

        # Try to handle NaNs for eventrate. These are handled differently from platform to platform (i.e. MSVC libraries are screwed). This platform-dependent fix is not needed in later versions of python. So, drop this in the near future!
        if str(self.eventrate) in ['-1.#IND', '1.#INF']:
            self.eventrate = 0

        # Only bits 0-19 are defined, zero the rest to make sure
        self.trigger_pattern &= 2 ** 20 - 1

        self.datetime = datetime.datetime(gps_year, gps_month, gps_day,
                                          gps_hour, gps_minute, gps_second)

        # Length of a single trace
        l = self.length / 2

        # Read out and save traces and calculated trace parameters
        self.mas_stdev1, self.mas_stdev2, self.mas_baseline1, \
        self.mas_baseline2, self.mas_npeaks1, self.mas_npeaks2, \
        self.mas_pulseheight1, self.mas_pulseheight2, self.mas_int1, \
        self.mas_int2, mas_tr1, mas_tr2 = \
        self.unpackSeqMessage('>8H2L%ds%ds' % (l, l))

        self.mas_tr1 = compress(self.unpack_trace(mas_tr1))
        self.mas_tr2 = compress(self.unpack_trace(mas_tr2))

        # Read out and save slave data as well, if available
        if self.num_devices > 1:
            self.slv_stdev1, self.slv_stdev2, self.slv_baseline1, \
            self.slv_baseline2, self.slv_npeaks1, self.slv_npeaks2, \
            self.slv_pulseheight1, self.slv_pulseheight2, self.slv_int1, \
            self.slv_int2, slv_tr1, slv_tr2 = \
            self.unpackSeqMessage('>8H2L%ds%ds' % (l, l))

            self.slv_tr1 = compress(self.unpack_trace(slv_tr1))
            self.slv_tr2 = compress(self.unpack_trace(slv_tr2))

    #--------------------------End of unpackMessage--------------------------#

    def unpack_trace(self, raw_trace):
        """    Unpack a trace
            Traces are stored in a funny way. We have a 12-bit ADC, so two
            datapoints can (and are) stored in 3 bytes. This function unravels
            traces again.

            DF: I'm wondering: does the LabVIEW program work hard to accomplish
            this? If so, why do we do this in the first place? The factor 1.5
            in storage space is hardly worth it, especially considering the
            fact that this is only used in the temporary buffer.

            DF: This is legacy code. I've never tried to understand it and will
            certainly not touch it until I do.
        """

        n = len(raw_trace)
        if n % 3 != 0:
            #return None
            raise Exception("Blob length is not divisible by 3!")
        a = struct.unpack("%dB" % (n), raw_trace)
        trace = []
        for i in xrange(0, n, 3):
            trace.append((a[i] << 4) + ((a[i + 1] & 240) >> 4))
            trace.append(((a[i + 1] & 15) << 8) + a[i + 2])
        trace_str = ""
        for i in trace:
            trace_str += str(i) + ","

        return trace_str

    #--------------------------End of unpack_trace--------------------------#
