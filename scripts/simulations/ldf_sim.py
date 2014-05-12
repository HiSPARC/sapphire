""" HiSPARC detector simulation

    This simulation takes an Extended Air Shower simulation ground
    particles file and uses that to simulate numerous showers hitting a
    HiSPARC detector station.  Only data of one shower is used, but by
    randomly selecting points on the ground as the position of a station,
    the effect of the same shower hitting various positions around the
    station is simulated.

"""
from __future__ import division

from math import pi

import tables

from sapphire import clusters
from sapphire.simulations import KascadeLdfSimulation


DATAFILE = 'data.h5'

N = 10000


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.open_file(DATAFILE, 'w')

    cluster = clusters.SingleStation()
    simulation = KascadeLdfSimulation(cluster, data, '/ldfsim/exact', R=60, N=N)
    simulation.run(max_theta=pi / 3)

    simulation = KascadeLdfSimulation(cluster, data, '/ldfsim/gauss_10', R=60, N=N, gauss=.1, trig_threshold=.9)
    simulation.run(max_theta=pi / 3)

    simulation = KascadeLdfSimulation(cluster, data, '/ldfsim/gauss_20', R=60, N=N, gauss=.2, trig_threshold=.8)
    simulation.run(max_theta=pi / 3)

    simulation = KascadeLdfSimulation(cluster, data, '/ldfsim/poisson', R=60, N=N, use_poisson=True)
    simulation.run(max_theta=pi / 3)

    simulation = KascadeLdfSimulation(cluster, data, '/ldfsim/poisson_gauss_20', R=60, N=N, use_poisson=True, gauss=.2, trig_threshold=.5)
    simulation.run(max_theta=pi / 3)
