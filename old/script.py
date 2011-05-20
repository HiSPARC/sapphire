import os, sys
sys.path.append(os.path.expanduser('~/work/HiSPARC/software/bzr/shower'))
import aires
import tables

from numpy import *

DATA_FILE = 'data-e15.h5'


class Particle(tables.IsDescription):
    id = tables.UInt32Col()
    pid = tables.Int8Col()
    core_distance = tables.Float32Col()
    polar_angle = tables.Float32Col()
    x = tables.Float32Col()
    y = tables.Float32Col()
    arrival_time = tables.Float32Col()
    energy = tables.Float32Col()


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

def store_simulation_data(data, group, file):
    sim = aires.SimulationData('sim', file)
    group = data.createGroup('/showers', group, 'Shower Data')
    gammas = data.createTable(group, 'gammas', Particle,
                              'Gammas')
    muons = data.createTable(group, 'muons', Particle,
                             'Muons and anti-muons')
    electrons = data.createTable(group, 'electrons', Particle,
                                 'Electrons and positrons')
    leptons = data.createTable(group, 'leptons', Particle,
                               'Electrons, positrons, muons and anti-muons')
    lepgammas = data.createTable(group, 'lepgammas', Particle,
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
    raise Exception("Warning!! Overwriting ALL data!!! (comment me)")
    data = tables.openFile(DATA_FILE, 'w')
    data.createGroup('/', 'showers', 'Simulated showers')
    store_simulation_data(data, 'zenith0', 'showere15.grdpcles')
    store_simulation_data(data, 'zenith23', 'showere15-angle.grdpcles')
    store_simulation_data(data, 'zenith5', 'showere15-angle-5.grdpcles')
    store_simulation_data(data, 'zenith35', 'showere15-angle-35.grdpcles')
    store_simulation_data(data, 'zenith40', 'showere15-angle-40.grdpcles')
    store_simulation_data(data, 'zenith45', 'showere15-angle-45.grdpcles')
    store_simulation_data(data, 'zenith60', 'showere15-angle-60.grdpcles')
    store_simulation_data(data, 'zenith80', 'showere15-angle-80.grdpcles')
