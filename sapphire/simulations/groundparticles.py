from math import pi, sin, cos, sqrt
import warnings

import numpy as np

from .base import BaseSimulation


class GroundParticlesSimulation(BaseSimulation):

    def __init__(self, corsikafile, max_core_distance, *args, **kwargs):
        super(GroundParticlesSimulation, self).__init__(*args, **kwargs)

        self.corsikafile = corsikafile
        self.max_core_distance = max_core_distance

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :returns: dictionary with shower parameters: core_pos
                  (x, y-tuple) and azimuth.

        """
        warnings.warn("Read some parameters from CORSIKA data!")

        R = self.max_core_distance

        # This is the fastest implementation I could come up with.  I
        # timed several permutations of numpy / math, and tried a Monte
        # Carlo method in which I pick x and y in a square (fast) and then
        # determine if they fall inside a circle or not (surprisingly
        # slow, because of an if-statement, and despite some optimizations
        # suggested by HM).

        for i in range(self.N):
            r = sqrt(np.random.uniform(0, R ** 2))
            phi = np.random.uniform(-pi, pi)
            alpha = np.random.uniform(-pi, pi)

            x = r * cos(phi)
            y = r * sin(phi)

            shower_parameters = {'core_pos': (x, y), 'azimuth': alpha}
            self._prepare_cluster_for_shower(shower_parameters)
            yield shower_parameters

    def _prepare_cluster_for_shower(self, shower_parameters):
        """Prepare the cluster object for the simulation of a shower.

        Rotate and translate the cluster so that (0, 0) coincides with the
        shower core position and 'East' coincides with the shower azimuth
        direction.

        """
        x, y = shower_parameters['core_pos']
        alpha = shower_parameters['azimuth']

        # rotate the core position around the original cluster center
        xp = x * cos(-alpha) - y * sin(-alpha)
        yp = x * sin(-alpha) + y * cos(-alpha)

        self.cluster.set_xyalpha_coordinates(-xp, -yp, -alpha)

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower."""

        # implement this!
        observables = None

        return observables

    def simulate_trigger(self, station_observables):
        """Simulate a trigger response."""

        return True
