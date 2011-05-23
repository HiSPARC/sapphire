""" HiSPARC detector simulation

    This simulation takes an Extended Air Shower simulation ground
    particles file and uses that to simulate numerous showers hitting a
    HiSPARC detector station.  Only data of one shower is used, but by
    randomly selecting points on the ground as the position of a station,
    the effect of the same shower hitting various positions around the
    station is simulated.

    For historical reasons, the detector simulation generates a CSV file
    as output.  This CSV file is then read and the data is stored in the
    HDF5 file.  Finally, the simulation results are analyzed to determine
    observables and these are also stored in the HDF5 file.

"""
from __future__ import division

import tables
import csv
import numpy as np
from math import pi, sqrt, sin, cos, atan2, ceil, isinf
import gzip
import progressbar as pb

from cluster_definition import cluster, DETECTOR_SIZE


DATAFILE = 'data-e15.h5'


class StationEvent(tables.IsDescription):
    """Store information about the station during an event"""
    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    r_x = tables.Float32Col()
    phi_y = tables.Float32Col()
    alpha = tables.Float32Col()

class ParticleEvent(tables.IsDescription):
    """Store information about the particles hitting a detector"""
    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    detector_id = tables.UInt8Col()
    pid = tables.Int8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    time = tables.Float32Col()
    energy = tables.Float32Col()

class ObservableEvent(tables.IsDescription):
    """Store information about the observables of an event"""
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
    d1 = tables.Float32Col()
    d2 = tables.Float32Col()
    d3 = tables.Float32Col()
    d4 = tables.Float32Col()


def generate_positions(R, num):
    """Generate positions and an orientation uniformly on a circle

    :param R: radius of circle
    :param num: number of positions to generate

    :return: r, phi, alpha

    """
    for i in range(num):
        phi, alpha = np.random.uniform(-pi, pi, 2)
        r = np.sqrt(np.random.uniform(0, R ** 2))
        yield r, phi, alpha

def get_station_coordinates(station, r, phi, alpha):
    """Calculate coordinates of a station given cluster coordinates

    :param station: station definition
    :param r, phi: polar coordinates of cluster center
    :param alpha: rotation of cluster

    :return: x, y, alpha; coordinates and rotation of station relative to
        absolute coordinate system

    """
    X = r * cos(phi)
    Y = r * sin(phi)

    sx, sy = station['position']
    xp = sx * cos(alpha) - sy * sin(alpha)
    yp = sx * sin(alpha) + sy * cos(alpha)

    x = X + xp
    y = Y + yp
    angle = alpha + station['angle']

    return x, y, angle

def get_station_particles(station, data, X, Y, alpha):
    """Return all particles hitting a station

    :param station: station definition
    :param data: HDF5 particle dataset
    :param X, Y: coordinates of station center
    :param alpha: rotation angle of station

    :return: list of detectors containing list of particles
    :rtype: list of lists

    """
    particles = []

    for detector in station['detectors']:
        x, y, orientation = detector
        particles.append(get_detector_particles(data, X, Y, x, y,
                                                orientation, alpha))
    return particles

def get_detector_particles(data, X, Y, x, y, orientation, alpha=None):
    """Return all particles hitting a single detector

    Given a HDF5 table containing information on all simulated particles
    and coordinates and orientation of a detector, search for all
    particles which have hit the detector.

    :param data: table containing particle data
    :param X, Y: X, Y coordinates of center of station
    :param x, y: x, y coordinates of center of detector relative to
        station center 
    :param orientation: either 'UD' or 'LR', for up-down or left-right
        detector orientations, relative to station
    :param alpha: rotation angle of entire station

    :return: list of particles which have hit the detector

    """
    c = get_detector_corners(X, Y, x, y, orientation, alpha)

    # determine equations describing detector edges
    b11, line1, b12 = get_line_boundary_eqs(c[0], c[1], c[2])
    b21, line2, b22 = get_line_boundary_eqs(c[1], c[2], c[3])

    # particles satisfying all equations are inside the detector
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

    :param p0, p1: (x, y) tuples on the same line
    :param p2: (x, y) tuple on the parallel line

    :return: value1, equation, value2, such that points satisfying value1
        < equation < value2 are between the parallel lines

    Example::

        >>> get_line_boundary_eqs((0, 0), (1, 1), (0, 2))
        (0.0, 'y - 1.000000 * x', 2.0)

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
    """Get the x, y coordinates of the detector corners

    :param X, Y: X, Y coordinates of center of station
    :param x, y: x, y coordinates of center of detector relative to
        station center 
    :param orientation: either 'UD' or 'LR', for up-down or left-right
        detector orientations, relative to station
    :param alpha: rotation angle of entire station

    :return: x, y coordinates of detector corners
    :rtype: list of (x, y) tuples

    """
    dx = DETECTOR_SIZE[0] / 2
    dy = DETECTOR_SIZE[1] / 2

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

    return [(X + x, Y + y) for x, y in corners]

def do_simulation(cluster, particles, data, dst, R, N):
    """Perform a simulation

    :param cluster: definition of all stations in the cluster
    :param particles: the HDF5 dataset containing the particles
    :param data: the HDF5 file
    :param dst: the HDF5 destination to store results
    :param R: maximum distance of shower to center of cluster
    :param N: number of simulations to perform

    """
    s_events = data.createTable(dst, 'stations', StationEvent)
    p_events = data.createTable(dst, 'particles', ParticleEvent)

    progress = pb.ProgressBar(maxval=N, widgets=[pb.Percentage(),
                                                 pb.Bar(), pb.ETA()])
    for event_id, (r, phi, alpha) in \
        progress(enumerate(generate_positions(R, N))):
        write_header(s_events, event_id, 0, r, phi, alpha)
        for station_id, station in enumerate(cluster):
            x, y, beta = get_station_coordinates(station, r, phi, alpha)
            write_header(s_events, event_id, station_id, x, y, beta)

            plist = get_station_particles(station, particles, x, y, beta)
            write_detector_particles(p_events, event_id, station_id,
                                     plist)

    s_events.flush()
    p_events.flush()

def write_header(table, event_id, station_id, r_x, phi_y, alpha):
    row = table.row
    row['id'] = event_id
    row['station_id'] = station_id
    row['r_x'] = r_x
    row['phi_y'] = phi_y
    row['alpha'] = alpha
    row.append()

def write_detector_particles(table, event_id, station_id, plist):
    row = table.row
    for detector_id, detector in enumerate(plist):
        for particle in detector:
            row['id'] = event_id
            row['station_id'] = station_id
            row['detector_id'] = detector_id
            row['pid'] = particle['pid']
            row['r'] = particle['core_distance']
            row['phi'] = particle['polar_angle']
            row['time'] = particle['arrival_time']
            row['energy'] = particle['energy']
            row.append()
    table.flush()

def analyze_results(hdffile, tablename, showertablename):
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

    showertable = data.getNode('/showers', showertablename)
    table = data.createTable('/analysis', tablename, StationEvent,
                             "Analyzed simulation data")
    sim = iter(data.getNode('/simulations', tablename))
    row = table.row

    dR = 3.
    rs = linspace(0 + dR, 100, 1000)
    dens = [len(showertable.readWhere('(%f <= core_distance) & '
                                      '(core_distance < %f)' % (R - dR,
                                                                R + dR))) /
            (pi * ((R + dR) ** 2 - (R - dR) ** 2)) for R in rs]

    try:
        event = sim.next()
        while True:
            assert event['id'] % 10 == 0
            row['id'] = event['id']
            row['r'] = event['r']
            R = event['r']
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
                    D = dens[rs.searchsorted(R) - 1]
                    row['d1'], row['d2'], row['d3'], row['d4'] = D, D, D, D
                    row.append()
                    break
                else:
                    t[idx - 1].append(event['time'])
    except StopIteration:
        row['t1'], row['t2'], row['t3'], row['t4'] = [min(x) if len(x)
                                                      else nan for x in t]
        row['n1'], row['n2'], row['n3'], row['n4'] = [len(x) for x in t]
        D = dens[rs.searchsorted(R) - 1]
        row['d1'], row['d2'], row['d3'], row['d4'] = D, D, D, D
        row.append()

    table.flush()
    data.close()


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'a')

    if '/simulations' not in data:
        data.createGroup('/', 'simulations', 'Detector Simulations')
    if 'zenith0' in data.root.simulations:
        data.removeNode('/simulations', 'zenith0', recursive=True)
    data.createGroup('/simulations', 'zenith0')

    particles = data.getNode('/showers/zenith0', 'leptons')
    dst = data.getNode('/simulations', 'zenith0')

    do_simulation(cluster, particles, data, dst, 100., N=1000)
