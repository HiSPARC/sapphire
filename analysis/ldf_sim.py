""" HiSPARC detector simulation

    This simulation takes an Extended Air Shower simulation ground
    particles file and uses that to simulate numerous showers hitting a
    HiSPARC detector station.  Only data of one shower is used, but by
    randomly selecting points on the ground as the position of a station,
    the effect of the same shower hitting various positions around the
    station is simulated.

"""
from __future__ import division

import tables
import os.path
import sys
import textwrap

import clusters
from simulations import KascadeLdfSimulation


DATAFILE = 'data.h5'


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'w')

    cluster = clusters.SimpleCluster()
    simulation = KascadeLdfSimulation(cluster, data, '/ldfsim', R=200,
                                      N=100000)
    simulation._Ne = 10 ** 6
    simulation.run()
