""" Simulation classes

    Define simulations as a class, so they can be subclassed.  This way,
    only bits of simulation logic need to be overridden.

"""
import tables
import os.path
import progressbar as pb
import textwrap
import numpy as np
from math import pi, sin, cos, atan2, sqrt, isinf

import storage


class BaseSimulation(object):

    """Base class

    This class is intended to be subclassed.  It provides no functionality, but
    shows the interface you can expect from simulations.

    """

    def __init__(self, cluster, data, output, R, N, force=False):
        """Simulation initialization

        :param cluster: BaseCluster (or derived) instance
        :param data: the HDF5 file
        :param output: name of the destination group to store results
        :param R: maximum distance of shower to center of cluster
        :param N: number of simulations to perform
        :param force: if True, ignore initialization errors, due to
            missing ground particles or existing previous simulations.
            Only use this if you really know what you're doing!

        """
        self.cluster = cluster
        self.data = data
        self.R = R
        self.N = N

        head, tail = os.path.split(output)
        try:
            self.output = self.data.createGroup(head, tail,
                                                createparents=True)
        except tables.NodeError:
            if force:
                self.output = self.data.getNode(head, tail)
            else:
                raise RuntimeError("Cancelling simulation; %s already exists?"
                                   % output)

    def run(self, positions=None):
        """Run a simulation

        This is the code which performs the simulation.  It creates a list
        of positions, creates all necessary tables and performs the
        simulation.

        :param positions: if given, use these coordinates instead of
            generating new ones

        """
        self._run_welcome_msg()

        if not positions:
            positions = self.generate_positions()

        # Example simulation loop
        #
        #for event_id, (r, phi, alpha) in enumerate(positions):
        #    self.simulate_event(event_id, r, phi, alpha)

        self._run_exit_msg()

    def _run_welcome_msg(self):
        """Print a welcome message at start of simulation run"""

        print 74 * '-'
        print textwrap.dedent("""Running simulation

                              Output destination: %s

                              Maximum core distance of cluster center:   %f m
                              Number of cluster positions in simulation: %d
                              """ % (self.output._v_pathname, self.R, self.N))

    def _run_exit_msg(self):
        print 74 * '-'
        print

    def write_observables(self, row, station, t):
        """Write observables from a single event

        :param row: PyTables row object
        :param station: station object
        :param t: 4-tuple of lists containing the arrival times of the
                  particles in the four detectors, e.g.
                  ([13.3, 8.6, 33.2], [], [3.4], [8.7, 3.6])

        """
        r = station['r']
        phi = station['phi']
        row['id'] = station['id']
        row['station_id'] = station['station_id']
        row['r'] = r
        row['phi'] = phi
        row['x'] = r * cos(phi)
        row['y'] = r * sin(phi)
        row['alpha'] = station['alpha']
        row['N'] = sum([1 if u else 0 for u in t])
        row['t1'], row['t2'], row['t3'], row['t4'] = \
            [min(u) if len(u) else -999 for u in t]
        row['n1'], row['n2'], row['n3'], row['n4'] = \
            [len(u) for u in t]
        row.append()

    def write_coincidence(self, row, event, N):
        """Write coincidence information

        :param row: PyTables row object
        :param event: simulated shower event information
        :param N: number of stations which triggered

        """
        r = event['r']
        phi = event['phi']
        row['id'] = event['id']
        row['N'] = N
        row['r'] = r
        row['phi'] = phi
        row['x'] = r * cos(phi)
        row['y'] = r * sin(phi)
        row.append()
