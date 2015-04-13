"""Perform simulations of CORSIKA air showers on a cluster of stations

This simulation uses a HDF5 file created from a CORSIKA simulation with
the `store_corsika_data` script. The shower is 'thrown' on the cluster
with random core positions and azimuth angles.

Example usage::

    import tables

    from sapphire.simulations.groundparticles import GroundParticlesSimulation
    from sapphire.clusters import ScienceParkCluster

    data = tables.open_file('/tmp/test_groundparticle_simulation.h5', 'w')
    cluster = ScienceParkCluster()

    sim = GroundParticlesSimulation('corsika.h5', 500, cluster, data, '/', 10)
    sim.run()

"""
from math import pi, sin, cos, sqrt

import numpy as np
import tables

from .detector import HiSPARCSimulation, ErrorlessSimulation
from ..utils import pbar


class GroundParticlesSimulation(HiSPARCSimulation):

    def __init__(self, corsikafile_path, max_core_distance, *args, **kwargs):
        """Simulation initialization

        :param corsikafile_path: path to the corsika.h5 file containing
                                 the groundparticles.
        :param max_core_distance: maximum distance of shower core to
                                  center of cluster.

        """
        super(GroundParticlesSimulation, self).__init__(*args, **kwargs)

        self.corsikafile = tables.open_file(corsikafile_path, 'r')
        self.groundparticles = self.corsikafile.get_node('/groundparticles')
        self.max_core_distance = max_core_distance

    def __del__(self):
        self.finish()

    def finish(self):
        """Clean-up after simulation"""

        self.corsikafile.close()

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :returns: dictionary with shower parameters: core_pos
                  (x, y-tuple) and azimuth.

        """
        r = self.max_core_distance
        giga = int(1e9)

        event_header = self.corsikafile.get_node_attr('/', 'event_header')
        event_end = self.corsikafile.get_node_attr('/', 'event_end')
        corsika_parameters = {'zenith': event_header.zenith,
                              'size': event_end.n_electrons_levels,
                              'energy': event_header.energy,
                              'particle': event_header.particle}

        for i in pbar(range(self.N)):
            x, y = self.generate_core_position(r)

            # Subtract Corsika shower azimuth from desired shower azimuth
            # make it fit in (-pi, pi] to get rotation angle of the cluster.
            shower_azimuth = self.generate_azimuth()
            alpha = shower_azimuth - event_header.azimuth
            if alpha > pi:
                alpha -= 2 * pi
            elif alpha <= -pi:
                alpha += 2 * pi

            shower_parameters = {'ext_timestamp': (giga + i) * giga,
                                 'core_pos': (x, y),
                                 'azimuth': shower_azimuth}
            shower_parameters.update(corsika_parameters)
            self._prepare_cluster_for_shower(x, y, alpha)

            yield shower_parameters

    def _prepare_cluster_for_shower(self, x, y, alpha):
        """Prepare the cluster object for the simulation of a shower.

        Rotate and translate the cluster so that (0, 0) coincides with the
        shower core position and 'East' coincides with the shower azimuth
        direction.

        :param x,y: position of shower core relative to cluster origin in m.
        :param alpha: angle the cluster needs to be rotated in radians.

        """
        # rotate the core position around the original cluster center
        xp = x * cos(-alpha) - y * sin(-alpha)
        yp = x * sin(-alpha) + y * cos(-alpha)

        self.cluster.set_coordinates(-xp, -yp, 0, -alpha)

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Checks if leptons have passed a detector. If so, it returns the number
        of leptons in the detector and the arrival time of the first lepton
        passing the detector.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         the observables will be determined.
        :param shower_parameters: dictionary with the shower parameters.

        """
        particles = self.get_particles_in_detector(detector)
        n_detected = len(particles)

        if n_detected:
            mips = self.simulate_detector_mips_for_particles(particles)
            particles['t'] += self.simulate_signal_transport_time(n_detected)
            first_signal = particles['t'].min() + detector.offset
            observables = {'n': mips,
                           't': self.simulate_adc_sampling(first_signal)}
        else:
            observables = {'n': 0., 't': -999}

        return observables

    def simulate_detector_mips_for_particles(self, particles):
        """Simulate the detector signal for particles

        :param particles: particle rows with the p_[x, y, z]
                          components of the particle momenta.

        """
        # determination of lepton angle of incidence
        theta = np.arccos(abs(particles['p_z']) /
                          np.sqrt(particles['p_x'] ** 2 +
                                  particles['p_y'] ** 2 +
                                  particles['p_z'] ** 2))
        n = len(particles)
        mips = self.simulate_detector_mips(n, theta)

        return mips

    def simulate_trigger(self, detector_observables):
        """Simulate a trigger response.

        :param detector_observables: dictionary containing the
                                     observables of one detector.
        :returns: True if at least 2 detectors detect at least one particle,
                  False otherwise.

        """
        detectors_hit = [True for observables in detector_observables
                         if observables['n'] > 0]

        if sum(detectors_hit) >= 2:
            return True
        else:
            return False

    def simulate_gps(self, station_observables, shower_parameters, station):
        """Simulate gps timestamp.

        :param station_observables: dictionary containing the observables
                                    of the station.
        :param shower_parameters: dictionary with the shower parameters.
        :param station: :class:`~sapphire.clusters.Station` for which
                         to simulate the gps timestamp.
        :returns: station_observables updated with gps timestamp and
                  trigger time.

        """
        arrival_times = [station_observables['t%d' % id]
                         for id in range(1, 5)
                         if station_observables.get('n%d' % id, -1) > 0]

        if len(arrival_times) > 1:
            trigger_time = sorted(arrival_times)[1]

            ext_timestamp = shower_parameters['ext_timestamp']
            ext_timestamp += int(trigger_time + station.gps_offset +
                                 self.simulate_gps_uncertainty())
            timestamp = int(ext_timestamp / int(1e9))
            nanoseconds = int(ext_timestamp % int(1e9))

            gps_timestamp = {'ext_timestamp': ext_timestamp,
                             'timestamp': timestamp,
                             'nanoseconds': nanoseconds,
                             't_trigger': trigger_time}
            station_observables.update(gps_timestamp)

        return station_observables

    def get_particles_in_detector(self, detector):
        """Get particles that hit a detector.

        Particle ids 2, 3, 5, 6 are electrons and muons,
        id 4 is no longer used (were neutrino's).

        The detector is approximated by a square with a surface of 0.5
        square meter which is *not* correctly rotated.  In fact, during
        the simulation, the rotation of the detector is undefined.  This
        is faster than a more thorough implementation.

        *Detector height is ignored!*

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.

        """
        x, y = detector.get_xy_coordinates()
        detector_boundary = sqrt(.5) / 2.
        query = ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
                 ' & (particle_id >= 2) & (particle_id <= 6)' %
                 (x - detector_boundary, x + detector_boundary,
                  y - detector_boundary, y + detector_boundary))
        return self.groundparticles.read_where(query)


class OptimizeQueryGroundParticlesSimulation(GroundParticlesSimulation):
    """
    Test pytables.read_where() in-kernel query performance versus
    combination of pytables.read_where()+numpy.compress()
    """
    def get_particles_in_detector(self, detector):
        """Get particles that hit a detector.

        Particle ids 2, 3, 5, 6 are electrons and muons,
        id 4 is no longer used (were neutrino's).

        The detector is approximated by a square with a surface of 0.5
        square meter which is *not* correctly rotated.  In fact, during
        the simulation, the rotation of the detector is undefined.  This
        is faster than a more thorough implementation.

        *Detector height is ignored!*

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.

        """
        x, y = detector.get_xy_coordinates()
        detector_boundary = sqrt(.5) / 2.

        """
        old:
        query = ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
                 ' & (particle_id >= 2) & (particle_id <= 6)' %
                 (x - detector_boundary, x + detector_boundary,
                  y - detector_boundary, y + detector_boundary))
        """
        pytables_X_query = ('(x >= %f) & (x <= %f)' %
                 (x - detector_boundary, x + detector_boundary))

        X_result = np.array(self.groundparticles.read_where(pytables_X_query))

        numpy_Y_query = ((X_result['y'] >= y - detector_boundary) & (X_result['y'] <= y + detector_boundary))
        lepton_query =  ((X_result['particle_id'] >= 2) & (X_result['particle_id'] <= 6))

        return X_result.compress(numpy_Y_query & lepton_query)

class OptimizeQuery_ParticlesOnly_GroundParticlesSimulation(GroundParticlesSimulation):
    """
    Test pytables.read_where() in-kernel query performance versus
    combination of pytables.read_where()+numpy.compress()
    """
    def get_particles_in_detector(self, detector):
        """Get particles that hit a detector.

        Particle ids 2, 3, 5, 6 are electrons and muons,
        id 4 is no longer used (were neutrino's).

        The detector is approximated by a square with a surface of 0.5
        square meter which is *not* correctly rotated.  In fact, during
        the simulation, the rotation of the detector is undefined.  This
        is faster than a more thorough implementation.

        *Detector height is ignored!*

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.

        """
        x, y = detector.get_xy_coordinates()
        detector_boundary = sqrt(.5) / 2.

        """
        old:
        query = ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
                 ' & (particle_id >= 2) & (particle_id <= 6)' %
                 (x - detector_boundary, x + detector_boundary,
                  y - detector_boundary, y + detector_boundary))
        """
        pytables_XY_query = ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)' %
                 (x - detector_boundary, x + detector_boundary,
                  y - detector_boundary, y + detector_boundary))

        XY = np.array(self.groundparticles.read_where(pytables_XY_query))

        lepton_query =  ((XY['particle_id'] >= 2) & (XY['particle_id'] <= 6))
        return XY.compress(lepton_query)

class DetectorBoundarySimulation(GroundParticlesSimulation):

    """ More accuratly simulate the detection area of the detectors.

    Take the orientation of the detectors into account and use the
    exact detector boundaries. This requires a slightly more complex
    query which is a bit slower.

    """

    def get_particles_in_detector(self, detector):
        """Simulate the detector detection area accurately.

        First particles are filtered to see which fall inside a
        non-rotated square box around the detector (i.e. sides of 1.2m).
        For the remaining particles a more accurate query is used to see
        which actually hit the detector. The advantage of using the
        square is that column indexes can be used, which may speed up
        queries.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.

        """
        detector_boundary = 0.6

        x, y = detector.get_xy_coordinates()
        c = detector.get_corners()

        b11, line1, b12 = self.get_line_boundary_eqs(c[0], c[1], c[2])
        b21, line2, b22 = self.get_line_boundary_eqs(c[1], c[2], c[3])
        query = ("(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f) & "
                 "(b11 < %s) & (%s < b12) & (b21 < %s) & (%s < b22) & "
                 "(particle_id >= 2) & (particle_id <= 6)" %
                 (x - detector_boundary, x + detector_boundary,
                  y - detector_boundary, y + detector_boundary,
                  line1, line1, line2, line2))

        return self.groundparticles.read_where(query)

    def get_line_boundary_eqs(self, p0, p1, p2):
        """Get line equations using three points

        Given three points, this function computes the equations for two
        parallel lines going through these points.  The first and second
        point are on the same line, whereas the third point is taken to
        be on a line which runs parallel to the first.  The return value
        is an equation and two boundaries which can be used to test if a
        point is between the two lines.

        :param p0, p1: (x, y) tuples on the same line
        :param p2: (x, y) tuple on the parallel line

        :return: value1, equation, value2, such that points satisfying
            value1 < equation < value2 are between the parallel lines

        Example::

            >>> get_line_boundary_eqs((0, 0), (1, 1), (0, 2))
            (0.0, 'y - 1.000000 * x', 2.0)

        """
        (x0, y0), (x1, y1), (x2, y2) = p0, p1, p2

        # Compute the general equation for the lines
        if not (x0 == x1):
            # First, compute the slope
            a = (y1 - y0) / (x1 - x0)

            # Calculate the y-intercepts of both lines
            b1 = y0 - a * x0
            b2 = y2 - a * x2

            line = "y - %f * x" % a
        else:
            # line is exactly vertical
            line = "x"
            b1, b2 = x0, x2

        # And order the y-intercepts
        if b1 > b2:
            b1, b2 = b2, b1

        return b1, line, b2


class ParticleCounterSimulation(GroundParticlesSimulation):

    """Do not simulate mips, just count the number of particles."""

    def simulate_detector_mips(self, n, theta):
        """A mip for a mip, count number of particles in a detector."""

        return n


class FixedCoreDistanceSimulation(GroundParticlesSimulation):

    """Shower core at a fixed core distance (from cluster origin).

    :param core_distance: distance of shower core to center of cluster.

    """

    def generate_core_position(cls, R):
        """Generate a random core position on a circle

        :param R: Core distance, in meters.
        :returns: Random x, y position on the circle with radius R.

        """
        phi = np.random.uniform(-pi, pi)
        x = R * cos(phi)
        y = R * sin(phi)
        return x, y


class GroundParticlesSimulationWithoutErrors(ErrorlessSimulation,
                                             GroundParticlesSimulation):

    """This simulation does not simulate errors/uncertainties

    This results in perfect timing (first particle through detector)
    and particle counting for the detectors.

    """

    pass
