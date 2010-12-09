""" Build datasets containing particle timing information

    This script will analyze the simulation data, make selections, and
    build datasets containing the timing information of pairs of
    particles.  The first particle in each pair will have a particular
    distance to the shower core and the second particle of the pair will
    be within a particular distance of the first one (typically
    inter-detector distances).

    These datasets are saved to a PyTables file.

"""

import tables
import time
from math import atan2, sqrt, pi
from numpy import linspace

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

    print "Performing analysis for %s and %s with %f <= R < %f" % (p1, p2,
                                                                   r0, r1)
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
    kwargs = dict(shower=shower, outfile=outfile, p1='leptons',
                  p2='leptons')

    bins = linspace(0, 80, 6)
    print 'Bin edges:', bins
    for r0, r1 in zip(bins[:-1], bins[1:]):
        make_particle_tables(r0=r0, r1=r1, dr0=0., dr1=1.,
                             table='g_%d_%d' % (r0, r1), **kwargs)

if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    try:
        outfile
    except NameError:
        outfile = tables.openFile(OUTFILE, 'w')
    
    do_analysis(data.root.showers.zenith0, outfile)
