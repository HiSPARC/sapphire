""" Store AIRES simulation data in HDF5 file

    This module reads the AIRES binary ground particles file and stores
    each particle individually in a HDF5 file, using PyTables.  This file
    can then be used as input for the detector simulation.

"""
import os, sys
sys.path.append(os.path.expanduser('~/work/HiSPARC/software/bzr/shower'))
import aires
import tables
import os.path

from sapphire.storage import ShowerParticle

from numpy import *

DATA_FILE = 'data.h5'


def save_particle(row, p, id):
    row['id'] = id
    row['pid'] = p.code
    row['core_distance'] = 10 ** p.core_distance
    row['polar_angle'] = p.polar_angle
    row['arrival_time'] = p.arrival_time
    row['energy'] = 10 ** p.energy
    row['x'] = 10 ** p.core_distance * cos(p.polar_angle)
    row['y'] = 10 ** p.core_distance * sin(p.polar_angle)
    row.append()


def store_aires_data(data, group_name, file):
    try:
        group = create_group(data, group_name)
    except tables.NodeError:
        print '%s already exists, doing nothing' % group_name
        return

    print "Storing AIRES data (%s) in %s" % (file, group)

    if not os.path.exists(file):
        raise RuntimeError("File %s does not exist" % file)
    else:
        sim = aires.SimulationData('', file)

    for shower_num, shower in enumerate(sim.showers()):
        shower_group = data.create_group(group, 'shower_%d' % shower_num)
        print shower_group

        leptons = data.create_table(shower_group, 'leptons', ShowerParticle,
                                   'Electrons, positrons, muons and anti-muons')

        leptons_row = leptons.row
        for id, p in enumerate(shower.particles()):
            if p.name == 'muon' or p.name == 'anti-muon' or \
                p.name == 'electron' or p.name == 'positron':
                save_particle(leptons_row, p, id)
        leptons.flush()

    sim.close()


def create_group(data, group_name):
    head, tail = os.path.split(group_name)
    group = data.create_group(head, tail, createparents=True)
    return group


if __name__ == '__main__':
    data = tables.open_file(DATA_FILE, 'a')
    store_aires_data(data, '/showers/E_1PeV/zenith_0',
                           '../aires/showere15-angle-0.grdpcles')
    store_aires_data(data, '/showers/E_100TeV/zenith_0',
                           '../aires/test.grdpcles')
