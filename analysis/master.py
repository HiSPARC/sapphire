import os
import warnings

import tables

import store_aires_data


class Master(object):
    def __init__(self, data_filename):
        if os.path.exists(data_filename):
            warnings.warn("%s already exists, some steps are skipped" % data_filename)
        self.data = tables.openFile(data_filename, 'a')

    def main(self):
        for angle in [0, 5, 22.5, 35]:
            self.store_1PeV_data_for_angle(angle)

    def store_1PeV_data_for_angle(self, angle):
        angle_string = str(angle).replace('.', '_')
        group_name = '/showers/E_1PeV/zenith_%s' % angle_string
        filename = '../aires/showere15-angle-%s.grdpcles' % angle_string

        store_aires_data.store_aires_data(self.data, group_name, filename)


if __name__ == '__main__':
    master = Master('master.h5')
    master.main()
