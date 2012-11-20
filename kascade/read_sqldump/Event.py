"""
    This module creates different types of Events that are specified by the subclasses
"""

__author__ = "thevinh"
__date__ = "$16-sep-2009"


class Event():
    # the instantiation operation
    def __init__(self):
        # init variables here if needed
        self.datetime = 0
        self.uploadCode = 0
        self.data = 0
        self.nanoseconds = 0
        self.export_values = 0

    #--------------------------End of __init__--------------------------#

    def getEventData(self):
        pass

    def parseMessage(self):
        pass
