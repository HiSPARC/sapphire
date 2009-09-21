""" Create tables

    Create a HiSPARC / KASCADE data structure in a new pytables file

"""
import tables
import os
from hisparc.containers import *

class Error(Exception):
    """Base class for Error exceptions"""
    pass

class FileExistsError(Error):
    """Raised when a file already exists"""
    def __init__(self, value='File already exists'):
        self.value = value

    def __str__(self):
        return self.value

def create_tables(filename):
    """Create HiSPARC and KASCADE tables in data file

    This function will create a new HDF5 file (if it doesn't exist already)
    and creates a basic structure useful for combined HiSPARC / KASCADE
    data analysis.

    Arguments:
    filename            the filename of the HDF5 file

    Return values:
    file                handle to the opened file

    """

    if os.path.exists(filename):
        raise FileExistsError

    print "Creating a new PyTables data file... ",
    data = tables.openFile(filename, 'a', 'HiSPARC / KASCADE data')
    hisparc = data.createGroup('/', 'hisparc', 'HiSPARC data')
    kascade = data.createGroup('/', 'kascade', 'KASCADE data')
    coincidences = data.createGroup('/', 'coincidences',
                                    'HiSPARC / KASCADE coincidences')
    data.createTable(hisparc, 'events', HisparcEvent, 'HiSPARC events',
                     expectedrows=50000)
    data.createVLArray(hisparc, 'traces', tables.VLStringAtom(),
                       'HiSPARC event traces', expectedsizeinMB=100)
    data.createTable(kascade, 'events', KascadeEvent, 'KASCADE events',
                     expectedrows=5000)
    data.createTable(coincidences, 'events', Coincidence,
                     'Coincidence events', expectedrows=5000)

    data.close()
    print "done."
