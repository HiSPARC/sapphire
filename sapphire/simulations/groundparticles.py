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
import random

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
            observables = {'n': round(mips, 3),
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

class GroundParticlesGammaSimulation(GroundParticlesSimulation):
    """ Implement digitisation of gamma photons """

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Checks if particles have passed a detector. If so, it returns the number
        of particles in the detector and the arrival time of the first particle
        passing the detector.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         the observables will be determined.
        :param shower_parameters: dictionary with the shower parameters.

        """

        leptons, gammas = self.get_particles_in_detector(detector)
        n_leptons = len(leptons)
        n_gammas = len(gammas)

        if not n_leptons+n_gammas:
            return {'n': 0, 't': -999}

        if n_leptons:
            mips_lepton = self.simulate_detector_mips_for_particles(leptons)
            leptons['t'] += self.simulate_signal_transport_time(n_leptons)
            first_signal = leptons['t'].min() + detector.offset
        else:
            mips_lepton = 0

        if n_gammas:
            mips_gamma = self.simulate_detector_mips_for_gammas(gammas)
            gammas['t'] += self.simulate_signal_transport_time(n_gammas)
            first_signal = gammas['t'].min() + detector.offset
        else:
            mips_gamma = 0

        return {'n': mips_lepton+mips_gamma, 't': self.simulate_adc_sampling(first_signal) }

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

        X_query = ('(x >= %f) & (x <= %f)' %
                 (x - detector_boundary, x + detector_boundary))

        X_selection = self.groundparticles.read_where(X_query)

        """
        use numpy.compress() to filter on y-coordinates
        this is faster than including the y-query in the pytables.read_where()
           for small (up to >100Mb) corsika hdf5 datafiles
        """
        Y_query = (X_selection['y'] >= (y-detector_boundary)) & (X_selection['y'] <= (y+detector_boundary))

        XY = X_selection.compress(Y_query)

        # use numpy.compress() to split leptons and gamma's
        leptons = (XY['particle_id'] >= 2) & (XY['particle_id'] <= 6)
        gammas = (XY['particle_id']  == 1)

        return XY.compress(leptons), XY.compress(gammas)

    def simulate_detector_mips_for_gammas(self, particles):
        """Simulate the detector signal for gammas

        :param particles: particle rows with the p_[x, y, z]
                          components of the particle momenta.

        """

        p_gamma = np.sqrt(particles['p_x'] ** 2 + particles['p_y'] ** 2 +
                particles['p_z'] ** 2)

        # determination of lepton angle of incidence
        theta = np.arccos(abs(particles['p_z']) /
                          p_gamma)

        n = len(particles)
        mips = self.simulate_detector_mips_gammas(n, p_gamma, theta)

        return mips

    @classmethod
    def simulate_detector_mips_gammas(cls, n, p, theta):
        """
        :param n: number of gammas.
        :param p: the momentum of the gammas as float or array
        :param theta: angle of incidence of the gammas, as float or array.
        """

        scintilator_depth = 2.0 # cm

        max_E = 4.0 # 2 MeV per cm * 2cm scintilator depth
        MIP = 3.38 # MeV

        L_rad = 42.52 # cm (radiation length in scinitlator)
        l_pair = 9./7. * 42.52 # mean free path for pair production W.R. Leo (1987) p.56

        # W.R. Leo (1987) p 54
        # E photon energy [MeV]
        # return compton edge [MeV]
        def _compton_edge(E):

            electron_rest_mass_MeV = .5109989 # MeV

            gamma = E / electron_rest_mass_MeV

            return (E * 2 * gamma / (1 + 2*gamma) )

        def _compton_energy_transfer(E):

            """
            from lio-project/photons/electron_energy_distribution.py

            The energy transfered from photon to electron = T(E) * compton_edge()
                the transfer function T(E) is calculated from dsigma/dT and
                represented as a lookup table.

            numpy.searchsorted() is a binarysearch to find the correct row in the lookup table
              => the polynomial coefficients corresponding to the given Energy (E)

            this is ugly code, but fast: 77.4 microsec per loop

            :param E: photon energy [MeV]
            :returns: energy transfered to electron [MeV]
            """

            Energy_table = np.array([\
                     0.100000 ,\
                     0.127427 ,\
                     0.162378 ,\
                     0.206914 ,\
                     0.263665 ,\
                     0.335982 ,\
                     0.428133 ,\
                     0.545559 ,\
                     0.695193 ,\
                     0.885867 ,\
                     1.128838 ,\
                     1.438450 ,\
                     1.832981 ,\
                     2.335721 ,\
                     2.976351 ,\
                     3.792690 ,\
                     4.832930 ,\
                     6.158482 ,\
                     7.847600 ,\
                     10.000000])

            transfer_function_table = [\
                    [-0.095663 , 0.998190 , 0.042602],\
                    [-0.104635 , 1.008633 , 0.040506],\
                    [-0.109294 , 1.015132 , 0.038090],\
                    [-0.107136 , 1.015115 , 0.035430],\
                    [-0.095551 , 1.005781 , 0.032673],\
                    [-0.072295 , 0.984547 , 0.030032],\
                    [-0.036015 , 0.949579 , 0.027764],\
                    [0.013328 , 0.900254 , 0.026124],\
                    [0.074324 , 0.837384 , 0.025313],\
                    [0.144294 , 0.763123 , 0.025434],\
                    [0.219738 , 0.680594 , 0.026476],\
                    [0.296884 , 0.593392 , 0.028324],\
                    [0.372186 , 0.505081 , 0.030787],\
                    [0.442678 , 0.418822 , 0.033635],\
                    [0.506156 , 0.337141 , 0.036638],\
                    [0.561214 , 0.261844 , 0.039593],\
                    [0.607175 , 0.194044 , 0.042338],\
                    [0.643960 , 0.134252 , 0.044755],\
                    [0.671940 , 0.082499 , 0.046776],\
                    [0.691794 , 0.038463 , 0.048368],
                    [0.691794 , 0.038463 , 0.048368]]  # extra item E > 10

            idx = Energy_table.searchsorted(E, side = 'left')
            p = np.poly1d(transfer_function_table[idx])
            return p(np.random.random())*_compton_edge(E)

        def _compton_mean_free_path(E):
            """
            Interaction probability based on Klein Nisihina cross section
            W.R. Leo (1987) p 54

            fit from lio-project/photons/mean_free_path.py

            :param E: photon energy [MeV]
            :returns: mean free path [cm]
            """
            if (E < 1.):
                return -0.806 * (E**2) + 8.46 * E + 6.36
            elif (E < 30.):
                return -0.109 * (E**2) + 6.00 * E + 8.53
            else:
                return 3.* E + 47. # [cm]

        #p [eV] and E [MeV]
        E = p / 1.e6

        mips = 0
        for energy, angle in zip(E,theta):

            costheta = cos(angle)

            """
            Calculate interaction point in units of scinitlator depth.
            depth > 1.0: no interaction
            """

            depth_compton = random.expovariate(1/_compton_mean_free_path(energy))
            depth_pair = random.expovariate(1/l_pair)

            if ((depth_pair > scintilator_depth/costheta) & \
                        (depth_compton > scintilator_depth/costheta)):
                # no interaction
                continue

            #  interaction in scintilator
            if (depth_compton < depth_pair):
                """
                compton scattering
                """
                maximum_energy_deposit_in_MIPS = ((scintilator_depth-depth_compton)/scintilator_depth)*max_E/MIP/costheta
                energy_deposit_in_MIPS = _compton_energy_transfer(energy)/MIP
                extra_mips = np.minimum(maximum_energy_deposit_in_MIPS, energy_deposit_in_MIPS)  # maximise energy transfer per photon to 1 MIP/cm * depth

                mips += extra_mips

            elif (energy > 1.022):
                """
                pair production: Two "electrons"
                """
                maximum_energy_deposit_in_MIPS = ((scintilator_depth-depth_pair)/scintilator_depth)*max_E/MIP/costheta
                energy_deposit_in_MIPS = (energy - 1.022) / MIP # 1.022 MeV used for creation of two particles
                extra_mips = np.minimum(maximum_energy_deposit_in_MIPS, energy_deposit_in_MIPS)

                mips += extra_mips

        return mips


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
