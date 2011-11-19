""" Simulation classes

    Define simulations as a class, so they can be subclassed.  This way,
    only bits of simulation logic need to be overridden.

"""
import tables
import os.path
import textwrap
from math import pi, sin, cos
import numpy as np

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
        :param force: if True, ignore pre-existing simulations; they will be
            overwritten!

        """
        self.cluster = cluster
        self.data = data
        self.R = R
        self.N = N

        if output in data and not force:
            raise RuntimeError("Cancelling simulation; %s already exists?"
                               % output)
        elif output in data:
            data.removeNode(output, recursive=True)

        head, tail = os.path.split(output)
        self.output = data.createGroup(head, tail, createparents=True)
        self.observables = self.data.createTable(self.output, 'observables',
                                                 storage.SimulationEventObservables)
        self.coincidences = self.data.createTable(self.output, 'coincidences',
                                                  storage.Coincidence)
        self.c_index = self.data.createVLArray(self.output, 'c_index',
                                               tables.UInt32Atom())

        self.output._v_attrs.cluster = cluster

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
        #for event_id, (r, phi) in enumerate(positions):
        #    self.simulate_event(event_id, r, phi)

        self._run_exit_msg()

    def _run_welcome_msg(self):
        """Print a welcome message at start of simulation run"""

        print 74 * '-'
        print textwrap.dedent("""\
            Running simulation

            Output destination: %s

            Maximum core distance of cluster center:   %f m
            Number of cluster positions in simulation: %d
            """ % (self.output._v_pathname, self.R, self.N))

    def _run_exit_msg(self):
        print 74 * '-'
        print

    def generate_positions(self):
        """Generate positions and an orientation uniformly on a circle

        :return: r, phi, alpha

        """
        for i in range(self.N):
            phi = np.random.uniform(-pi, pi)
            r = np.sqrt(np.random.uniform(0, self.R ** 2))
            yield r, phi

    def write_observables(self, station, t):
        """Write observables from a single event

        :param station: station object
        :param t: 4-tuple of lists containing the arrival times of the
                  particles in the four detectors, e.g.
                  ([13.3, 8.6, 33.2], [], [3.4], [8.7, 3.6])

        """
        row = self.observables.row

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

    def write_coincidence(self, event, N):
        """Write coincidence information

        :param event: simulated shower event information
        :param N: number of stations which triggered

        """
        row = self.coincidences.row

        r = event['r']
        phi = event['phi']
        row['id'] = event['id']
        row['N'] = N
        row['r'] = r
        row['phi'] = phi
        row['x'] = r * cos(phi)
        row['y'] = r * sin(phi)
        row.append()
