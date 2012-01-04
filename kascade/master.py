""" Read KASCADE data into analysis file

    This module reads a KASCADE ascii datafile and stores the contents
    inside the pytables analysis file

"""
import tables
import logging
from hisparc import kascade
from hisparc.containers import KascadeEvent


class Master(object):
    def __init__(self, data_filename, kascade_filename):
        self.data = tables.openFile(data_filename, 'a')
        self.kascade_filename = kascade_filename

    def main(self):
        self.read_and_store_kascade_data()

    def read_and_store_kascade_data(self):
        """Read KASCADE data into analysis file"""

        if 'kascade' in self.data.root:
            logging.info("KASCADE event group already exists, skipping "
                         "reconstruction")
            return
        else:
            self.data.createGroup('/', 'kascade', "KASCADE data")
            self.data.createTable('/kascade', 'events', KascadeEvent,
                                  "KASCADE events")
            kascade.helper(self.data.root.hisparc.cluster_kascade.station_601.events,
                           self.data.root.kascade.events,
                           self.kascade_filename)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    master = Master('kascade.h5', 'HiSparc-new.dat.gz')
    master.main()
