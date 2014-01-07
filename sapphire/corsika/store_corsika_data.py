#!/usr/bin/env python

""" Store CORSIKA simulation data in HDF5 file

    This module reads the CORSIKA binary ground particles file and stores
    each particle individually in a HDF5 file, using PyTables.  This file
    can then be used as input for the detector simulation.

"""
import os
import sys

import tables
from numpy import *

from sapphire import corsika

SOURCE_FILE = '/Users/niekschultheiss/corsika/corsika-74000/run/DAT000002'
DEST_FILE = '/Users/niekschultheiss/data/data.h5'


class GroundParticles(tables.IsDescription):
    particle_name = tables.StringCol(30)
    r = tables.Float32Col()
    phi = tables.Float32Col()
    x = tables.Float32Col()
    y = tables.Float32Col()
    t = tables.Float32Col()
    p_x = tables.Float32Col()
    p_y = tables.Float32Col()
    p_z = tables.Float32Col()


def save_particle(row, p):
    row['particle_name'] = p.particle
    row['r'] = p.r
    row['phi'] = p.phi
    row['x'] = p.x
    row['y'] = p.y
    row['t'] = p.t
    row['p_x'] = p.p_x
    row['p_y'] = p.p_y
    row['p_z'] = p.p_z
    row.append()


def store_corsika_data(source, group_name, destination):
    try:
        table = destination.createTable('/', group_name, GroundParticles,
                                 'All groundparticles')
    except tables.NodeError:
        print '%s already exists, doing nothing' % group_name
        return

    print "Storing CORSIKA data (%s) in %s" % (source._filename,
                                               destination.filename)

    particle_row = table.row

    for event in source.get_events():
        for particle in event.get_particles():
            save_particle(particle_row, particle)

    table.flush()


def main():
    # Source
    corsika_data = corsika.CorsikaFile(SOURCE_FILE)
    corsika_data.check()

    # Destination      
    hdf_data = tables.openFile(DEST_FILE, 'a')

    store_corsika_data(corsika_data, 'groundparticles', hdf_data)

    hdf_data.close()


if __name__ == '__main__':
    main()
