""" Read KASCADE data into analysis file

    This module reads a KASCADE ascii datafile and stores the contents
    inside the pytables analysis file

"""
import tables
from hisparc import kascade
from hisparc.containers import KascadeEvent


DATAFILE = 'kascade.h5'
KASCADEFILE = 'HiSparc-new.dat.gz'


def read_kascade_data():
    """Read KASCADE data into analysis file"""

    data = tables.openFile(DATAFILE, 'a')

    if 'kascade' in data.root:
        raise RuntimeError("KASCADE event group already exists!")
    else:
        data.createGroup('/', 'kascade', "KASCADE data")
        data.createTable('/kascade', 'events', KascadeEvent,
                         "KASCADE events")
        kascade.helper(data.root.hisparc.cluster_kascade.station_601.events,
                       data.root.kascade.events, KASCADEFILE)


if __name__ == '__main__':
    read_kascade_data()
