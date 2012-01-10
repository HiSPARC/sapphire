import logging
from math import pi

import tables

from sapphire.kascade import StoreKascadeData, KascadeCoincidences
from sapphire.storage import KascadeEvent
from sapphire.analysis.process_events import ProcessIndexedEvents
from sapphire import clusters
from sapphire.analysis.direction_reconstruction import KascadeDirectionReconstruction


class Master(object):
    hisparc_group = '/hisparc/cluster_kascade/station_601'
    kascade_group = '/kascade'

    def __init__(self, data_filename, kascade_filename):
        self.data = tables.openFile(data_filename, 'a')
        self.kascade_filename = kascade_filename

    def main(self):
        self.store_cluster_instance()
        self.read_and_store_kascade_data()
        self.search_for_coincidences()
        self.process_events()
        self.reconstruct_direction()

    def store_cluster_instance(self):
        group = self.data.getNode(self.hisparc_group)

        if 'cluster' not in group._v_attrs:
            cluster = clusters.SingleStation()
            cluster.set_xyalpha_coordinates(65., 20.82, pi)

            group._v_attrs.cluster = cluster

    def read_and_store_kascade_data(self):
        """Read KASCADE data into analysis file"""

        try:
            kascade = StoreKascadeData(self.data, self.hisparc_group,
                                       self.kascade_group,
                                       self.kascade_filename)
        except RuntimeError, msg:
            print msg
            return
        else:
            kascade.read_and_store_data()

    def search_for_coincidences(self):
        hisparc = self.hisparc_group
        kascade = self.kascade_group

        try:
            coincidences = KascadeCoincidences(self.data, hisparc, kascade)
        except RuntimeError, msg:
            print msg
            return
        else:
            print "Searching for coincidences"
            coincidences.search_coincidences(timeshift= -13.180220188,
                                             dtlimit=1e-3)
            print "Storing coincidences"
            coincidences.store_coincidences()
            print "Done."

    def process_events(self):
        c_index = self.data.getNode(self.kascade_group, 'c_index')
        index = c_index[:, 1]

        process = ProcessIndexedEvents(self.data, self.hisparc_group,
                                       index)
        try:
            process.process_and_store_results()
        except RuntimeError, msg:
            print msg
            return

    def reconstruct_direction(self):
        reconstruction = KascadeDirectionReconstruction(self.data, '/reconstructions', overwrite=True)
        reconstruction.reconstruct_angles(self.hisparc_group, self.kascade_group)


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    master = Master('kascade.h5', 'HiSparc-new.dat.gz')
    master.main()
