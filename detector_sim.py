""" HiSPARC detector simulation

    This simulation takes an Extended Air Shower simulation ground
    particles file and uses that to simulate numerous showers hitting a
    HiSPARC detector station.  Only data of one shower is used, but by
    randomly selecting points on the ground as the position of a station,
    the effect of the same shower hitting various positions around the
    station is simulated.

    To enhance the statistics, (far) more positions in lower particle
    density regions are selected.

"""
import tables
import csv
import numpy as np
from numpy import pi, sqrt, arctan2, sin, cos

from pylab import *

#OUTFILE = 'simulation-e15.csv'
OUTFILE = 'tmp.csv'
DATAFILE = 'data-e15.h5'
HDFFILE = 'simulation-e15.h5'

D = .01
MIN = 3

RINGS = [(0, 4, 20, False), (4, 20, 10, False), (20, 40, 16, False),
         (40, 80, 30, False), (80, 80, 30, True)]
DETECTORS = [(0., 5.77, 'UD'), (0., 0., 'UD'), (-5., -2.89, 'LR'),
             (5., -2.89, 'LR')]
DETECTOR_SIZE = (-.25, .25, -.5, .5)


class ParticleEvent(tables.IsDescription):
    id = tables.UInt32Col()
    pid = tables.Int8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    time = tables.Float32Col()
    energy = tables.Float32Col()

class StationEvent(tables.IsDescription):
    id = tables.UInt32Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    times = tables.Float32Col(shape=4)


def simulate_positions(r0, r1=1., density=1., iscorner=False):
    """Simulate station positions on a ring with uniform density"""

    if not iscorner:
        num = ceil(density * pi * (r1 ** 2 - r0 ** 2))
        if num < MIN:
            num = MIN
        print r0, r1, density, num
        return random_ring(r0, r1, num)
    else:
        num = ceil(density * ((2 * r0) ** 2 - pi * r0 ** 2))
        if num < MIN:
            num = MIN
        print r0, "Corner", density, num
        return random_corner(r0, num)

def random_ring(r0, r1, num):
    """Simulate positions uniformly on a ring"""

    phi = np.random.uniform(-pi, pi, num)
    r = sqrt(np.random.uniform(r0 ** 2, r1 ** 2, num))

    return r, phi

def random_corner(R, num):
    """Simulate positions in a square with an inscribed circle left out"""

    r_list, phi_list = [], []
    while len(r_list) < num:
        x, y = np.random.uniform(-R, R, 2)
        r = sqrt(x ** 2 + y ** 2)
        if r < R:
            continue
        else:
            phi = arctan2(y, x)
            r_list.append(r)
            phi_list.append(phi)
    return np.array(r_list), np.array(phi_list)

def plot_positions_test():
    """Test of the position generator"""

    N = 0
    figure()

    for r0, r1, density, iscorner in RINGS:
        r, phi = simulate_positions(r0, r1, density * D, iscorner)
        N += len(r)
        plot(r * sin(phi), r * cos(phi), '.', ms=1)

    axis('equal')
    title("Posities in shower")
    xlabel("Afstand (m)")
    ylabel("Afstand (m)")
    legend(["0 - 4 m", "4 - 20 m", "20 - 40 m", "40 - 80 m", "> 80 m"])

    print "Total:", N

def get_station_particles(data, r, phi):
    """Return all particles hitting a station"""

    X = r * cos(phi)
    Y = r * sin(phi)
    particles = []

    for detector in DETECTORS:
        x, y, orientation = detector
        particles.append(get_detector_particles(data, X + x, Y + y,
                                                orientation))
    return particles

def get_detector_particles(data, x, y, orientation):
    """Return all particles hitting a single detector"""

    x0, x1, y0, y1 = get_detector_corners(x, y, orientation)
    return data.readWhere("(x0 < x) & (x < x1) & (y0 < y) & (y < y1)")

def get_detector_corners(x, y, orientation):
    """Get the x, y coordinates of the detector corners"""

    dx0, dx1, dy0, dy1 = DETECTOR_SIZE

    if orientation == 'UD':
        return x + dx0, x + dx1, y + dy0, y + dy1
    elif orientation == 'LR':
        return x + dy0, x + dy1, y + dx0, y + dx1
    else:
        raise Exception("Unknown detector orientation: %s" % orientation)

def plot_detectors_test(data, r, phi):
    """Test the detector and particle positions"""

    figure()
    X = r * cos(phi)
    Y = r * sin(phi)

    for x, y, orient in DETECTORS:
        plot(X + x, Y + y, '.', ms=5.)
        x0, x1, y0, y1 = get_detector_corners(X + x, Y + y, orient)
        print x0, x1, y0, y1
        plot([x0, x0, x1, x1, x0], [y0, y1, y1, y0, y0])

    title("HiSPARC detectors")
    xlabel("Afstand (m)")
    ylabel("Afstand (m)")
    axis('equal')

    for p in get_station_particles(data, r, phi):
        plot(p[:]['x'], p[:]['y'], '.', ms=1.)

def detector_test(data):
    """Test the detector positions and measured particles"""

    plot_detectors_test(group, 0, 0)
    plot_detectors_test(group, 5, 0)
    plot_detectors_test(group, 5, .5 * pi)
    plot_detectors_test(group, 50, .25 * pi)

def do_simulation(data, density):
    """Perform a full simulation"""

    with open(OUTFILE, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        N = 0
        for r0, r1, rel_density, iscorner in RINGS:
            r_list, phi_list = simulate_positions(r0, r1,
                                                  rel_density * density,
                                                  iscorner)
            for event_id, (r, phi) in enumerate(zip(r_list, phi_list), N):
                save_event_header(writer, event_id, r, phi)
                particles = get_station_particles(data, r, phi)
                for scint_id, p in enumerate(particles, 1):
                    save_detector_particles(writer, event_id, scint_id, p)

            N = event_id + 1

def save_event_header(writer, event_id, r, phi):
    writer.writerow([event_id * 10, 0, r, phi, 0, 0])

def save_detector_particles(writer, event_id, scint_id, particles):
    for p in particles:
        writer.writerow([event_id * 10 + scint_id, p['pid'],
                         p['core_distance'], p['polar_angle'],
                         p['arrival_time'], p['energy']])

def store_results_in_tables(csvfile, hdffile):
    """Read a csv file with simulation data and store in PyTables"""

    data = tables.openFile(hdffile, 'w')
    data.createTable('/', 'sim', ParticleEvent,
                     "Detector simulation data")
    table = data.root.sim
    row = table.row

    with open(csvfile, 'r') as file:
        reader = csv.reader(file, delimiter='\t')
        for line in reader:
            (row['id'], row['pid'], row['r'], row['phi'], row['time'],
             row['energy']) = line
            row.append()
    table.flush()
    data.close()

def analyze_results(hdffile):
    """Analyze simulation results and deduce detector pulse times"""

    data = tables.openFile(hdffile, 'a')

    try:
        data.removeNode('/', 'analysis')
    except tables.NoSuchNodeError:
        pass

    data.createTable('/', 'analysis', StationEvent,
                     "Analyzed detector simulation data")

    sim = iter(data.root.sim)
    table = data.root.analysis
    row = table.row

    try:
        event = sim.next()
        while True:
            assert event['id'] % 10 == 0
            row['id'] = event['id']
            row['r'] = event['r']
            row['phi'] = event['phi']
            t = [[], [], [], []]
            while True:
                event = sim.next()
                idx = event['id'] % 10
                if idx == 0:
                    row['times'] = [min(x) if len(x) else nan for x in t]
                    row.append()
                    break
                else:
                    t[idx - 1].append(event['time'])
    except StopIteration:
        row['times'] = [min(x) if len(x) else nan for x in t]
        row.append()

    table.flush()
    data.close()


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')
    
    try:
        data2
    except NameError:
        data2 = tables.openFile(HDFFILE, 'r')

    group = data.root.showers.s1.leptons

    #plot_positions_test()
    #detector_test(group)
    #_ip.magic("time do_simulation(group, .0001)")
    #_ip.magic("time do_simulation(group, .0002)")
    #_ip.magic("time do_simulation(group, .0004)")
    #_ip.magic("time do_simulation(group, D)")

    do_simulation(group, .001)
    #store_results_in_tables(OUTFILE, HDFFILE)
    #analyze_results(HDFFILE)
