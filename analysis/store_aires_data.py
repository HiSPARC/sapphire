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

from storage import ShowerParticle

from numpy import *

DATA_FILE = 'data-e15.h5'


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

def store_aires_data(data, group, file):
    head, tail = os.path.split(group)
    try:
        data.createGroup(head, tail, createparents=True)
    except tables.NodeError:
        print "Ignoring", group, file
        print "(%s already exists?)" % group
        return

    print "Storing AIRES data (%s) in %s" % (file, group)

    sim = aires.SimulationData('', file)
    gammas = data.createTable(group, 'gammas', ShowerParticle,
                              'Gammas')
    muons = data.createTable(group, 'muons', ShowerParticle,
                             'Muons and anti-muons')
    electrons = data.createTable(group, 'electrons', ShowerParticle,
                                 'Electrons and positrons')
    leptons = data.createTable(group, 'leptons', ShowerParticle,
                               'Electrons, positrons, muons and anti-muons')
    lepgammas = data.createTable(group, 'lepgammas', ShowerParticle,
                                 'Electrons, positrons, muons, '
                                 'anti-muons and gammas')

    shower = [x for x in sim.showers()][0]

    muons_row = muons.row
    gammas_row = gammas.row
    electrons_row = electrons.row
    leptons_row = leptons.row
    lepgammas_row = lepgammas.row

    for id, p in enumerate(shower.particles()):
        if p.name == 'gamma':
            save_particle(gammas_row, p, id)
            save_particle(lepgammas_row, p, id)
        elif p.name == 'muon' or p.name == 'anti-muon':
            save_particle(muons_row, p, id)
            save_particle(leptons_row, p, id)
            save_particle(lepgammas_row, p, id)
        elif p.name == 'electron' or p.name == 'positron':
            save_particle(electrons_row, p, id)
            save_particle(leptons_row, p, id)
            save_particle(lepgammas_row, p, id)
    muons.flush()
    gammas.flush()
    electrons.flush()
    leptons.flush()
    lepgammas.flush()


if __name__ == '__main__':
    data = tables.openFile(DATA_FILE, 'a')
    store_aires_data(data, '/showers/E_1PeV/zenith_0',
                     '../aires/showere15-angle-0.grdpcles')
