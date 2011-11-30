import os
import warnings
import re

import tables

import store_aires_data
from simulations import GroundParticlesSimulation, QSubSimulation
import clusters


R = 100
N = 100000
N_CORES = 32


class Master(object):
    def __init__(self, data_filename):
        if os.path.exists(data_filename):
            warnings.warn("%s already exists, some steps are skipped" % data_filename)
        self.data = tables.openFile(data_filename, 'a')

    def main(self):
        self.store_shower_data()
        self.do_station_simulations()
        self.do_custom_simulations()

    def store_shower_data(self):
        for angle in [0, 5, 22.5, 35]:
            self.store_1PeV_data_for_angle(angle)

    def store_1PeV_data_for_angle(self, angle):
        angle_string = str(angle).replace('.', '_')
        group_name = '/showers/E_1PeV/zenith_%s' % angle_string
        filename = '../aires/showere15-angle-%s.grdpcles' % angle_string

        store_aires_data.store_aires_data(self.data, group_name, filename)

    def do_station_simulations(self):
        cluster = clusters.SingleStation()
        for angle in self.get_shower_angles_from_shower_data():
            for shower in self.get_showers_in_group(angle):
                self.perform_simulation(cluster, shower)

    def do_custom_simulations(self):
        angle = '/showers/E_1PeV/zenith_22_5'
        for station_size in [5, 20]:
            cluster = clusters.SingleStation(station_size)
            for shower in self.get_showers_in_group(angle):
                output_path = self.make_output_path_for_station_size_simulation(shower, station_size)
                self.perform_simulation(cluster, shower, output_path)

    def make_output_path_for_station_size_simulation(self, shower, station_size):
        path = shower._v_pathname.replace('/showers/', '/simulations/')
        return re.sub('(/zenith_[0-9_]+)/', r'\1_size%d/' % station_size, path)

    def perform_simulation(self, cluster, shower, output_path=None):
        if not output_path:
            output_path = shower._v_pathname.replace('/showers/', '/simulations/')
        shower_path = shower.leptons._v_pathname

        args = [cluster, self.data, shower_path, output_path]
        kwargs = {'R': R, 'N': N}

        if self.is_qsub_available():
            Simulation = QSubSimulation
            kwargs['N_cores'] = N_CORES
        else:
            Simulation = GroundParticlesSimulation

        try:
            sim = Simulation(*args, **kwargs)
        except RuntimeError, msg:
            print msg
            return
        else:
            sim.run()

    def get_shower_angles_from_shower_data(self):
        return self.data.listNodes('/showers/E_1PeV')

    def get_showers_in_group(self, group):
        return self.data.listNodes(group)

    def is_qsub_available(self):
        return os.path.exists('/usr/bin/qsub')


if __name__ == '__main__':
    master = Master('master.h5')
    master.main()
