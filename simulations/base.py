""" Simulation classes

    Define simulations as a class, so they can be subclassed.  This way,
    only bits of simulation logic need to be overridden.

"""
import tables
import os.path
import progressbar as pb
import numpy as np
from numpy import nan
from math import pi, sin, cos, atan2, sqrt, isinf

import storage


class BaseSimulation(object):

    """Base simulation class

    This class defines a detector simulation, taking a shower simulation
    as input.  This class can be overridden to, for example, parallelize
    the simulation.
    """

    def __init__(self, cluster, data, grdpcles, output, R, N):
        """Simulation initialization

        :param cluster: BaseCluster (or derived) instance
        :param data: the HDF5 file
        :param grdpcles: name of the dataset containing the ground particles
        :param output: name of the destination group to store results
        :param R: maximum distance of shower to center of cluster
        :param N: number of simulations to perform

        """
        self.cluster = cluster
        self.data = data
        self.R = R
        self.N = N

        try:
            self.grdpcles = data.getNode('/', grdpcles)
        except tables.NoSuchNodeError:
            raise RuntimeError("Cancelling simulation; %s not found in "
                               "tree." % grdpcles)

        head, tail = os.path.split(output)
        try:
            self.output = self.data.createGroup(head, tail,
                                                createparents=True)
        except tables.NodeError:
            raise RuntimeError("Cancelling simulation; %s already exists?"
                               % output)

    def generate_positions(self):
        """Generate positions and an orientation uniformly on a circle

        :return: r, phi, alpha

        """
        for i in range(self.N):
            phi, alpha = np.random.uniform(-pi, pi, 2)
            r = np.sqrt(np.random.uniform(0, self.R ** 2))
            yield r, phi, alpha

    def get_station_coordinates(self, station, r, phi, alpha):
        """Calculate coordinates of a station given cluster coordinates

        :param station: station definition
        :param r, phi: polar coordinates of cluster center
        :param alpha: rotation of cluster

        :return: x, y, alpha; coordinates and rotation of station relative to
            absolute coordinate system

        """
        X = r * cos(phi)
        Y = r * sin(phi)

        sx, sy = station.position
        xp = sx * cos(alpha) - sy * sin(alpha)
        yp = sx * sin(alpha) + sy * cos(alpha)

        x = X + xp
        y = Y + yp
        angle = alpha + station.angle

        return x, y, angle

    def get_station_particles(self, station, X, Y, alpha):
        """Return all particles hitting a station

        :param station: station definition
        :param X, Y: coordinates of station center
        :param alpha: rotation angle of station

        :return: list of detectors containing list of particles
        :rtype: list of lists

        """
        particles = []

        size = station.detector_size
        for detector in station.detectors:
            x, y, orientation = detector
            particles.append(self.get_detector_particles(X, Y, x, y, size,
                                                         orientation,
                                                         alpha))
        return particles

    def get_detector_particles(self, X, Y, x, y, size, orientation,
                               alpha):
        """Return all particles hitting a single detector

        Given a HDF5 table containing information on all simulated particles
        and coordinates and orientation of a detector, search for all
        particles which have hit the detector.

        :param X, Y: X, Y coordinates of center of station
        :param x, y: x, y coordinates of center of detector relative to
            station center 
        :param size: tuple (width, length) giving detector size
        :param orientation: either 'UD' or 'LR', for up-down or left-right
            detector orientations, relative to station
        :param alpha: rotation angle of entire station

        :return: list of particles which have hit the detector

        """
        c = self.get_detector_corners(X, Y, x, y, size, orientation,
                                      alpha)

        # determine equations describing detector edges
        b11, line1, b12 = self.get_line_boundary_eqs(c[0], c[1], c[2])
        b21, line2, b22 = self.get_line_boundary_eqs(c[1], c[2], c[3])

        # particles satisfying all equations are inside the detector
        return self.grdpcles.readWhere("(b11 < %s) & (%s < b12) & "
                                       "(b21 < %s) & (%s < b22)" % \
                                       (line1, line1, line2, line2))

    def get_line_boundary_eqs(self, p0, p1, p2):
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

    def get_detector_corners(self, X, Y, x, y, size, orientation, alpha):
        """Get the x, y coordinates of the detector corners

        :param X, Y: X, Y coordinates of center of station
        :param x, y: x, y coordinates of center of detector relative to
            station center 
        :param size: tuple (width, length) giving detector size
        :param orientation: either 'UD' or 'LR', for up-down or left-right
            detector orientations, relative to station
        :param alpha: rotation angle of entire station

        :return: x, y coordinates of detector corners
        :rtype: list of (x, y) tuples

        """
        dx = size[0] / 2
        dy = size[1] / 2

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

    def run(self):
        """Run a simulation

        This method is just an entry function.  It can easily be
        overridden without the need to rewrite parts of the simulation.

        In this form, it just generates positions and calls
        :meth:`_do_run`.

        """
        positions = self.generate_positions()
        self._do_run(positions)

    def _do_run(self, positions):
        """Perform the actual simulation

        This is the actual code which performs the simulation.  It takes a
        list or a generator of positions, creates all necessary tables and
        performs the simulation.

        :param positions: list or generator of the positions to be
            simulated

        """
        print 74 * '-'
        print """Running simulation

Ground particles:   %s
Output destination: %s

Maximum core distance of cluster center:   %f m
Number of cluster positions in simulation: %d
        """ % (self.grdpcles._v_pathname, self.output._v_pathname, self.R,
               self.N)

        headers = self.data.createTable(self.output, 'headers',
                                         storage.SimulationHeader)
        particles = self.data.createTable(self.output, 'particles',
                                         storage.ParticleEvent)

        progress = pb.ProgressBar(maxval=self.N, widgets=[pb.Percentage(),
                                                          pb.Bar(),
                                                          pb.ETA()])
        for event_id, (r, phi, alpha) in progress(enumerate(positions)):
            self.write_header(headers, event_id, 0, r, phi, alpha)
            for station_id, station in enumerate(self.cluster.stations, 1):
                x, y, beta = self.get_station_coordinates(station, r, phi,
                                                          alpha)
                # calculate station r, phi just to save it in header
                s_r = sqrt(x ** 2 + y ** 2)
                s_phi = atan2(y, x)
                self.write_header(headers, event_id, station_id, s_r,
                                  s_phi, beta)

                plist = self.get_station_particles(station, x, y, beta)
                self.write_detector_particles(particles, event_id,
                                              station_id, plist)

        headers.flush()
        particles.flush()

        print 74 * '-'
        print

    def write_header(self, table, event_id, station_id, r, phi, alpha):
        """Write simulation event header information to file

        :param table: HDF5 table
        :param event_id: simulation event id
        :param station_id: station id inside cluster, 0 for cluster header
        :param r, phi: r, phi for cluster or station position, both as
            absolute coordinates
        :param alpha: cluster rotation angle or station rotation angle

        """
        row = table.row
        row['id'] = event_id
        row['station_id'] = station_id
        row['r'] = r
        row['phi'] = phi
        row['alpha'] = alpha
        row.append()
        table.flush()

    def write_detector_particles(self, table, event_id, station_id, plist):
        """Write particles to file

        :param table: HDF5 table
        :param event_id: simulation event id
        :param station_id: station id inside cluster
        :param plist: list of detectors, containing list of particles

        """
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

    def store_observables(self):
        """Analyze simulation results and store derived data

        Loop through simulation result tables and find observables like the
        number of particles which hit detectors, as well as the arrival time
        of the first particle to hit a detector.  The number of detectors
        which have particles are also recorded.  Finally the per shower
        results from all stations are combined and stored as a coincidence
        event.

        To speed things up, pointers into the two simulation result tables are
        advanced row by row.  Event and station id's are continually checked
        to now when to break.  The flow is a bit complicated, but it is fast.

        """
        obs = self.data.createTable(self.output, 'observables',
                                    storage.ObservableEvent)
        coinc = self.data.createTable(self.output, 'coincidences',
                                      storage.CoincidenceEvent)

        headers = self.data.getNode(self.output, 'headers')
        particles = self.data.getNode(self.output, 'particles')

        print "Storing observables from %s" % self.output._v_pathname

        obs_row = obs.row
        coinc_row = coinc.row
        progress = pb.ProgressBar(maxval=len(headers),
                                  widgets=[pb.Percentage(), pb.Bar(),
                                           pb.ETA()]).start()
        headers = iter(headers)
        particles = iter(particles)

        # start with first rows initialized
        header = headers.next()
        particle = particles.next()
        # loop over events
        while True:
            assert header['station_id'] == 0
            # freeze header row for later use
            event = header.fetch_all_fields()

            # N = number of stations which trigger
            N = 0
            try:
                # loop over headers
                while True:
                    header = headers.next()
                    # if station_id == 0, we have a new event, not a station
                    if header['station_id'] == 0:
                        break
                    else:
                        t = [[], [], [], []]
                        # loop over particles
                        while True:
                            # check if particle matches current event/station
                            if particle['id'] != event['id'] or \
                                particle['station_id'] != \
                                header['station_id']:
                                break
                            else:
                                t[particle['detector_id']].append(particle['time'])
                                try:
                                    particle = particles.next()
                                except StopIteration:
                                    # Ran out of particles. Forcing invalid
                                    # id, so that higher-level loops can
                                    # finish business, but no new particles
                                    # will be processed
                                    particle['id'] = -1
                        self.write_observables(obs_row, header, t)
                        # trigger if Ndet hit >= 2
                        if sum([1 if u else 0 for u in t]) >= 2:
                            N += 1
            # StopIteration when we run out of headers
            except StopIteration:
                break
            finally:
                self.write_coincidence(coinc_row, event, N)
                progress.update(header.nrow + 1)
        progress.finish()

        obs.flush()
        coinc.flush()
        print

    def write_observables(self, row, station, t):
        row['id'] = station['id']
        row['station_id'] = station['station_id']
        row['r'] = station['r']
        row['phi'] = station['phi']
        row['alpha'] = station['alpha']
        row['N'] = sum([1 if u else 0 for u in t])
        row['t1'], row['t2'], row['t3'], row['t4'] = \
            [min(u) if len(u) else nan for u in t]
        row['n1'], row['n2'], row['n3'], row['n4'] = \
            [len(u) for u in t]
        row.append()

    def write_coincidence(self, row, event, N):
        row['id'] = event['id']
        row['N'] = N
        row['r'] = event['r']
        row['phi'] = event['phi']
        row['alpha'] = event['alpha']
        row.append()
