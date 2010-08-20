""" Build datasets containing particle timing information

    This script will analyze the simulation data, make selections, and
    build datasets containing the timing information of pairs of
    particles.  The first particle in each pair will have a particular
    distance to the shower core (typically 150 m) and the second particle
    of the pair will be within a particular distance of the first one
    (typically inter-detector distances).

    These datasets are saved to a PyTables file as well as a CSV file.

"""

import tables
import csv
import time
from math import atan2, sqrt, pi

DATAFILE = 'data-e15.h5'
OUTFILE = 'analysis-e15.h5'

class ParticlePair(tables.IsDescription):
    PID = tables.Int8Col()
    R = tables.Float32Col()
    PHI = tables.Float32Col()
    T = tables.Float32Col()
    E = tables.Float32Col()
    pid = tables.Int8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    t = tables.Float32Col()
    e = tables.Float32Col()


def make_particle_tables(shower, outfile, table, p1, p2, r0, r1, dr0, dr1,
                         q1='', q2=''):
    if q1:
        q1 = ' & ' + q1
    if q2:
        q2 = ' & ' + q2

    try:
        outfile.removeNode('/', table)
    except tables.NoSuchNodeError:
        pass

    table = outfile.createTable('/', table, ParticlePair,
                                "Analysis of %s and %s with core "
                                "distance from %.1f to %.1f meters and a "
                                "surrounding circle between %.1f and "
                                "%.1f meters" % (p1, p2, r0, r1, dr0, dr1))

    row = table.row

    dr0sq = dr0 ** 2
    dr1sq = dr1 ** 2

    print "Performing analysis for %s and %s" % (p1, p2)
    t0 = time.time()
    N, M = 0, 0
    for p in getattr(shower, p1).where('(r0 <= core_distance) & '
                                       '(core_distance < r1)' + q1):
        x0, y0 = p['x'], p['y']
        for q in getattr(shower, p2).where(
                '(dr0sq <= ((x - x0) ** 2 + (y - y0) ** 2)) & '
                '(((x - x0) ** 2 + (y - y0) ** 2) < dr1sq)' + q2):
            if p['id'] != q['id']:
                x, y = q['x'], q['y']

                phi = atan2(y - y0, x - x0) - p['polar_angle']
                if phi > pi:
                    phi -= 2 * pi
                elif phi < -pi:
                    phi += 2 * pi

                row['PID'] = p['pid']
                row['R'] = p['core_distance']
                row['PHI'] = p['polar_angle']
                row['T'] = p['arrival_time']
                row['E'] = p['energy']
                row['pid'] = q['pid']
                row['r'] = sqrt((x - x0) ** 2 + (y - y0) ** 2)
                row['phi'] = phi
                row['t'] = q['arrival_time']
                row['e'] = q['energy']
                row.append()
                M += 1
            else:
                N += 1
    table.flush()
    t1 = time.time()
    print "Rejected: %d self-referencing pairs" % N
    print "Analysis took %.1f seconds (%d pairs)" % (t1 - t0, M)

def do_analysis(shower, outfile):
    kwargs = dict(shower=shower, outfile=outfile, p1='lepgammas',
                  p2='lepgammas')

    make_particle_tables(r0=50., r1=60., dr0=0., dr1=1.,
                         table='g50_60_0', **kwargs)
    make_particle_tables(r0=50., r1=60., dr0=5.2, dr1=6.2,
                         table='g50_60_6', **kwargs)
    make_particle_tables(r0=50., r1=60., dr0=9.5, dr1=10.5,
                         table='g50_60_10', **kwargs)
    
    make_particle_tables(r0=20., r1=25., dr0=0., dr1=1.,
                         table='g20_25_0', **kwargs)
    make_particle_tables(r0=20., r1=25., dr0=5.2, dr1=6.2,
                         table='g20_25_6', **kwargs)
    make_particle_tables(r0=20., r1=25., dr0=9.5, dr1=10.5,
                         table='g20_25_10', **kwargs)

    make_particle_tables(r0=10., r1=15., dr0=0., dr1=1.,
                         table='g10_15_0', **kwargs)
    make_particle_tables(r0=10., r1=15., dr0=5.2, dr1=6.2,
                         table='g10_15_6', **kwargs)
    make_particle_tables(r0=10., r1=15., dr0=9.5, dr1=10.5,
                         table='g10_15_10', **kwargs)

if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    try:
        outfile
    except NameError:
        outfile = tables.openFile(OUTFILE, 'a')
    
    shower = data.root.showers.s1

    do_analysis(shower, outfile)
