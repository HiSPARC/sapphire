from math import pi, sin, cos, sqrt
import warnings

import numpy as np

from .base import BaseSimulation


class GroundParticlesSimulation(BaseSimulation):

    def __init__(self, corsikafile, max_core_distance, *args, **kwargs):
        super(GroundParticlesSimulation, self).__init__(*args, **kwargs)

        self.corsikafile = corsikafile
        self.groundparticles = corsikafile.getNode('/groundparticles')
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
        """Simulate detector response to a shower.

        Checks if leptons have passed a detector. If so, it returns the number
        of leptons in the detector and the arrival time of the first lepton
        passing the detector.

        The detector is approximated by a square with a surface of 0.5
        square meter which is *not* correctly rotated.  In fact, during
        the simulation, the rotation of the detector is undefined.  This
        might be faster than a more thorough implementation.

        """
        detector_boundary = 0.3535534
        x, y = detector.get_xy_coordinates()

        # particle ids 2, 3, 5, 6 are electrons and muons, and id 4 is no
        # longer used (were neutrino's).
        query = ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
                 ' & (particle_id >= 2) & (particle_id <= 6)' %
                 (x - detector_boundary, x + detector_boundary,
                  y - detector_boundary, x + detector_boundary))
        detected = [row['t'] for row in self.groundparticles.where(query)]

        observables = {'n': len(detected), 't': min(detected)}

        return observables

    def simulate_trigger(self, station_observables):
        """Simulate a trigger response.

        :returns: True if at least 2 detectors detect at least one particle,
                  False otherwise.

        """
        trigger = False
        hitted_plates = 0

        for observables in station_observables:
            if observables['n'] > 0:
                hitted_plates += 1

        if hitted_plates > 1:
            trigger = True

        return trigger
