import tables
import logging
from sapphire.kascade import StoreKascadeData, KascadeCoincidences
from sapphire.storage import KascadeEvent
from sapphire.analysis.process_events import ProcessEvents


class Master(object):
    hisparc_group = '/hisparc/cluster_kascade/station_601'

    def __init__(self, data_filename, kascade_filename):
        self.data = tables.openFile(data_filename, 'a')
        self.kascade_filename = kascade_filename

    def main(self):
        self.read_and_store_kascade_data()
        self.search_for_coincidences()
        self.process_events()

    def read_and_store_kascade_data(self):
        """Read KASCADE data into analysis file"""

        try:
            kascade = StoreKascadeData(self.data, self.hisparc_group,
                                       '/kascade', self.kascade_filename)
        except RuntimeError, msg:
            print msg
            return
        else:
            kascade.read_and_store_data()

    def search_for_coincidences(self):
        hisparc = self.hisparc_group
        kascade = '/kascade'

        try:
            coincidences = KascadeCoincidences(self.data, hisparc, kascade)
        except RuntimeError, msg:
            print msg
            return
        else:
            coincidences.search_coincidences(timeshift= -13.180220188, dtlimit=1e-3)
            coincidences.store_coincidences()

    def process_events(self):
        process = ProcessEvents(self.data, self.hisparc_group,
                                overwrite=True)
        process.process_and_store_results()


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)

    master = Master('kascade.h5', 'HiSparc-new.dat.gz')
    master.main()
