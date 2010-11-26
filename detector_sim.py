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
from math import pi, sqrt, sin, cos, atan2, ceil, isinf
import gzip

from multiprocessing import Process

from pylab import *

DATAFILE = 'data-e15.h5'

D = 1.
MIN = 3

RINGS = [(0, 4, 20, False), (4, 20, 10, False), (20, 40, 16, False),
         (40, 80, 30, False), (80, 80, 30, True)]

DETECTOR_SIZE = (.25, .5)


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
    alpha = tables.Float32Col()
    t1 = tables.Float32Col()
    t2 = tables.Float32Col()
    t3 = tables.Float32Col()
    t4 = tables.Float32Col()
    n1 = tables.UInt16Col()
    n2 = tables.UInt16Col()
    n3 = tables.UInt16Col()
    n4 = tables.UInt16Col()


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
    r = np.sqrt(np.random.uniform(r0 ** 2, r1 ** 2, num))

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
            phi = atan2(y, x)
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

def get_station_particles(data, r, phi, alpha=None):
    """Return all particles hitting a station"""

    X = r * cos(phi)
    Y = r * sin(phi)
    particles = []

    for detector in DETECTORS:
        x, y, orientation = detector
        particles.append(get_detector_particles(data, X, Y, x, y,
                                                orientation, alpha))
    return particles

def get_detector_particles(data, X, Y, x, y, orientation, alpha=None):
    """Return all particles hitting a single detector"""

    c = get_detector_corners(X, Y, x, y, orientation, alpha)

    b11, line1, b12 = get_line_boundary_eqs(c[0], c[1], c[2])
    b21, line2, b22 = get_line_boundary_eqs(c[1], c[2], c[3])

    return data.readWhere("(b11 < %s) & (%s < b12) & "
                          "(b21 < %s) & (%s < b22)" % (line1, line1,
                                                       line2, line2))

def get_line_boundary_eqs(p0, p1, p2):
    """Get line equations using three points

    Given three points, this function computes the equations for two
    parallel lines going through these points.  The first and second point
    are on the same line, whereas the third point is taken to be on a
    line which runs parallel to the first.  The return value is an
    equation and two boundaries which can be used to test if a point is
    between the two lines.

    """
    (x0, y0), (x1, y1), (x2, y2) = p0, p1, p2

    # First, compute the slope
    a = (y1 - y0) / (x1 - x0)

    # Calculate the y-intercepts of both lines
    b1 = y0 - a * x0
    b2 = y2 - a * x2

    # Compute the general equation for the lines
    if not isinf(a):
        line = "y - %f * x" % a
    else:
        # line is exactly vertical
        line = "x"
        b1, b2 = x0, x2

    # And order the y-intercepts
    if b1 > b2:
        b1, b2 = b2, b1

    return b1, line, b2

def get_detector_corners(X, Y, x, y, orientation, alpha=None):
    """Get the x, y coordinates of the detector corners"""

    dx, dy = DETECTOR_SIZE

    if orientation == 'UD':
        corners = [(x - dx, y - dy), (x + dx, y - dy), (x + dx, y + dy),
                   (x - dx, y + dy)]
    elif orientation == 'LR':
        corners = [(x - dy, y - dx), (x + dy, y - dx), (x + dy, y + dx),
                   (x - dy, y + dx)]
    else:
        raise Exception("Unknown detector orientation: %s" % orientation)

    if alpha is not None:
        sina = sin(alpha)
        cosa = cos(alpha)
        corners = [[x * cosa - y * sina, x * sina + y * cosa] for x, y in
                   corners]

    return [[X + x, Y + y] for x, y in corners]

def plot_detectors_test(data, r, phi, alpha=None):
    """Test the detector and particle positions"""

    figure()
    X = r * cos(phi)
    Y = r * sin(phi)

    for x, y, orient in DETECTORS:
        plot(X + x, Y + y, '.', ms=5.)
        corners = get_detector_corners(X, Y, x, y, orient, alpha)
        plot([u for u, v in corners] + [corners[0][0]],
             [v for u, v in corners] + [corners[0][1]])

    title("HiSPARC detectors")
    xlabel("Afstand (m)")
    ylabel("Afstand (m)")
    axis('equal')

    for p in get_station_particles(data, r, phi, alpha):
        plot(p[:]['x'], p[:]['y'], '.', ms=1.)

def detector_test(data):
    """Test the detector positions and measured particles"""

    for theta in linspace(-pi, pi, 10):
        plot_detectors_test(group, 0, 0, theta)
        plot_detectors_test(group, 5, 0, theta)
        #plot_detectors_test(group, 5, .5 * pi, theta)
        #plot_detectors_test(group, 50, .25 * pi, theta)

def do_simulation(data, stationsize, outfile, density, use_alpha=0.):
    """Perform a full simulation"""

    dataf = tables.openFile(DATAFILE, 'r')
    data = dataf.getNode('/showers', data)

    global DETECTORS
    DETECTOR_SPACING = stationsize
    l = DETECTOR_SPACING / 2.
    x = l / 3. * sqrt(3)
    DETECTORS = [(0., 2 * x, 'UD'), (0., 0., 'UD'), (-l, -x, 'LR'),
                 (l, -x, 'LR')]

    with gzip.open(outfile, 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        N = 0
        for r0, r1, rel_density, iscorner in RINGS:
            r_list, phi_list = simulate_positions(r0, r1,
                                                  rel_density * density,
                                                  iscorner)
            if use_alpha is True:
                alpha_list = np.random.uniform(-pi, pi, len(r_list))
            else:
                alpha_list = len(r_list) * [use_alpha]

            for event_id, (r, phi, alpha) in enumerate(zip(r_list,
                                                           phi_list,
                                                           alpha_list),
                                                       N):
                save_event_header(writer, event_id, r, phi, alpha)
                particles = get_station_particles(data, r, phi, alpha)
                for scint_id, p in enumerate(particles, 1):
                    save_detector_particles(writer, event_id, scint_id, p)

            N = event_id + 1
    dataf.close()

def save_event_header(writer, event_id, r, phi, alpha):
    # ID + 0, 0, r, phi, alpha, 0
    writer.writerow([event_id * 10, 0, r, phi, alpha, 0])

def save_detector_particles(writer, event_id, scint_id, particles):
    # ID + scintnum, pid, r, phi, t, E
    for p in particles:
        writer.writerow([event_id * 10 + scint_id, p['pid'],
                         p['core_distance'], p['polar_angle'],
                         p['arrival_time'], p['energy']])

def store_results_in_tables(csvfile, hdffile):
    """Read a csv file with simulation data and store in PyTables"""

    data = tables.openFile(hdffile, 'a')
    try:
        data.createGroup('/', 'simulations', "Detector simulation data")
    except tables.NodeError:
        pass

    tablename = csvfile.replace('detsim-', '').replace(
                    '.csv.gz', '').replace('-', '_')

    try:
        data.removeNode('/simulations', tablename)
    except tables.NoSuchNodeError:
        pass

    table = data.createTable('/simulations', tablename, ParticleEvent,
                             "Detector simulation data")
    row = table.row

    file = gzip.open(csvfile, 'r')
    reader = csv.reader(file, delimiter='\t')
    for line in reader:
        (row['id'], row['pid'], row['r'], row['phi'], row['time'],
         row['energy']) = line
        row.append()

    file.close()
    table.flush()
    data.close()

def analyze_results(hdffile, tablename):
    """Analyze simulation results and deduce detector pulse times"""

    data = tables.openFile(hdffile, 'a')
    try:
        data.createGroup('/', 'analysis', "Analyzed simulation data")
    except tables.NodeError:
        pass

    try:
        data.removeNode('/analysis', tablename)
    except tables.NoSuchNodeError:
        pass

    table = data.createTable('/analysis', tablename, StationEvent,
                             "Analyzed simulation data")
    sim = iter(data.getNode('/simulations', tablename))
    row = table.row

    try:
        event = sim.next()
        while True:
            assert event['id'] % 10 == 0
            row['id'] = event['id']
            row['r'] = event['r']
            row['phi'] = event['phi']
            # header has alpha in place of time
            row['alpha'] = event['time']
            t = [[], [], [], []]
            while True:
                event = sim.next()
                idx = event['id'] % 10
                if idx == 0:
                    row['t1'], row['t2'], row['t3'], row['t4'] = \
                        [min(x) if len(x) else nan for x in t]
                    row['n1'], row['n2'], row['n3'], row['n4'] = \
                        [len(x) for x in t]
                    row.append()
                    break
                else:
                    t[idx - 1].append(event['time'])
    except StopIteration:
        row['t1'], row['t2'], row['t3'], row['t4'] = [min(x) if len(x)
                                                      else nan for x in t]
        row['n1'], row['n2'], row['n3'], row['n4'] = [len(x) for x in t]
        row.append()

    table.flush()
    data.close()

def do_full_simulation():
    DENSITY = 1.
    ps = []
    for group, stationsize, outfile in [
        ('zenith0/leptons', 10., 'detsim-angle-0.csv.gz'),
        ('zenith5/leptons', 10., 'detsim-angle-5.csv.gz'),
        ('zenith23/leptons', 5., 'detsim-angle-23-size5.csv.gz'),
        ('zenith23/leptons', 10., 'detsim-angle-23.csv.gz'),
        ('zenith23/leptons', 20., 'detsim-angle-23-size20.csv.gz'),
        ('zenith35/leptons', 10., 'detsim-angle-35.csv.gz'),
        ('zenith80/leptons', 10., 'detsim-angle-80.csv.gz')]:

        p = Process(target=do_simulation, kwargs=dict(data=group,
                                                      stationsize=stationsize,
                                                      outfile=outfile,
                                                      density=DENSITY,
                                                      use_alpha=True))
        p.start()
        ps.append(p)
    for p in ps:
        p.join()

def store_full_results():
    store_results_in_tables('detsim-angle-0.csv.gz', DATAFILE)
    store_results_in_tables('detsim-angle-5.csv.gz', DATAFILE)
    store_results_in_tables('detsim-angle-23.csv.gz', DATAFILE)
    store_results_in_tables('detsim-angle-23-size5.csv.gz', DATAFILE)
    store_results_in_tables('detsim-angle-23-size20.csv.gz', DATAFILE)
    store_results_in_tables('detsim-angle-35.csv.gz', DATAFILE)
    store_results_in_tables('detsim-angle-80.csv.gz', DATAFILE)

def analyze_full_results():
    analyze_results(DATAFILE, 'angle_0')
    analyze_results(DATAFILE, 'angle_5')
    analyze_results(DATAFILE, 'angle_23')
    analyze_results(DATAFILE, 'angle_23_size5')
    analyze_results(DATAFILE, 'angle_23_size20')
    analyze_results(DATAFILE, 'angle_35')
    analyze_results(DATAFILE, 'angle_80')


if __name__ == '__main__':
    do_full_simulation()
    store_full_results()
    analyze_full_results()
