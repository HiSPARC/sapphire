from __future__ import division

import os
import tables

import clusters
from simulations import GroundParticlesSimulation, QSubSimulation


DATAFILE = 'data-e15.h5'


if __name__ == '__main__':
    data = tables.openFile(DATAFILE, 'a')

    if '/simulations' in data:
        print "Removing previous simulations."
        data.removeNode('/simulations', recursive=True)

    sim = 'E_1PeV/zenith_0'
    cluster = clusters.SimpleCluster()
    simulation = GroundParticlesSimulation(cluster, data,
                                           os.path.join('/showers', sim,
                                                        'leptons'),
                                           os.path.join('/simulations', sim),
                                           R=100, N=10)
    simulation.run()
