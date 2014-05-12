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
from simulations import GroundParticlesSimulation, QSubSimulation


DATAFILE = 'data.h5'


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.open_file(DATAFILE, 'a')

    if '/simulations' in data:
        print
        print textwrap.dedent("""\
            WARNING: previous simulations exist and will be overwritten
            Continue? (answer 'yes'; anything else will exit)""")
        try:
            inp = raw_input()
        except KeyboardInterrupt:
            inp = 'Ctrl-C'

        if inp.lower() == 'yes':
            data.remove_node('/simulations', recursive=True)
        else:
            print
            print "Aborting!"
            sys.exit(1)

    sim = 'E_1PeV/zenith_0'
    cluster = clusters.SimpleCluster()
    simulation = GroundParticlesSimulation(cluster, data,
                                           os.path.join('/showers', sim,
                                                        'leptons'),
                                           os.path.join('/simulations',
                                                        sim),
                                           R=100, N=100)
    simulation.run()
