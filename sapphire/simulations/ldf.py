"""
Lateral distribution functions that can be used for simulating particle
densities and for fitting to data.

Example usage::

    import tables

    from sapphire.simulations.ldf import BaseLdfSimulation
    from sapphire.clusters import ScienceParkCluster

    data = tables.open_file('/tmp/test_ldf_simulation.h5', 'w')
    cluster = ScienceParkCluster()

    sim = NkgLdfSimulation(max_core_distance=400, cluster=cluster,
                           datafile=data, N=200)
    sim.run()

"""
from scipy.special import gamma
from numpy import pi, sin, cos, sqrt, random

from .detector import HiSPARCSimulation, ErrorlessSimulation
from ..utils import pbar


class BaseLdfSimulation(HiSPARCSimulation):

    def __init__(self, max_core_distance, *args, **kwargs):
        """Simulation initialization

        :param max_core_distance: maximum distance of shower core to
                                  center of cluster.

        """
        super(BaseLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = BaseLdf()
        self.max_core_distance = max_core_distance

        # The cluster is not moved, so detector positions can be stored.
        for station in self.cluster.stations:
            for detector in station.detectors:
                detector.xy_coordinates = detector.get_xy_coordinates()

    def generate_shower_parameters(self):
        """Generate shower parameters, i.e. core position

        For the simple LDF only the core position is relevant. It
        assumes the shower to come from the Zenith.

        :returns: dictionary with shower parameters: core_pos
                  (x, y-tuple).

        """
        r = self.max_core_distance
        giga = int(1e9)

        for i in pbar(range(self.N)):
            shower_parameters = {'ext_timestamp': (giga + i) * giga,
                                 'azimuth': self.generate_azimuth(),
                                 'zenith': 0.,
                                 'core_pos': self.generate_core_position(r),
                                 'size': None,
                                 'energy': None}

            yield shower_parameters

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower

        Get the mips in a detector from the LDF.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         the observables will be determined.
        :param shower_parameters: dictionary with the shower parameters.

        """
        n_detected = self.get_num_particles_in_detector(detector,
                                                        shower_parameters)
        theta = shower_parameters['zenith']

        if n_detected:
            mips = self.simulate_detector_mips(n_detected, theta)
            observables = {'n': mips}
        else:
            observables = {'n': 0.}
        return observables

    def get_num_particles_in_detector(self, detector, shower_parameters):
        """Get the number of particles in a detector

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         the number of particles will be determined.
        :param shower_parameters: dictionary with the shower parameters.

        """
        x, y = detector.xy_coordinates
        core_x, core_y = shower_parameters['core_pos']
        zenith = shower_parameters['zenith']
        azimuth = shower_parameters['azimuth']

        r = self.ldf.calculate_core_distance_from_coordinates_and_direction(
                x, y, core_x, core_y, zenith, azimuth)
        p_shower = self.ldf.calculate_ldf_value(r)
        p_ground = p_shower * cos(zenith)
        num_particles = random.poisson(p_ground * detector.get_area())

        return num_particles


class NkgLdfSimulation(BaseLdfSimulation):

    """Same as the BaseLdfSimulation but uses the NkgLdf as LDF"""

    def __init__(self, *args, **kwargs):
        super(NkgLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = NkgLdf()


class KascadeLdfSimulation(BaseLdfSimulation):

    """Same as the BaseLdfSimulation but uses the KascadeLdf as LDF"""

    def __init__(self, *args, **kwargs):
        super(KascadeLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = KascadeLdf()


class BaseLdf(object):

    def calculate_ldf_value(self, r):
        return 0.

    def calculate_core_distance_from_coordinates_and_direction(self,
                                                               x, y, x0, y0,
                                                               theta, phi):
        """Calculate core distance

        The core distance is the distance of the detector to the shower core,
        measured *on the shower front*.  For derivations, see logbook.

        :param x,y: detector position in m.
        :param x0,y0: shower core position in m.
        :param theta,phi: shower axis direction in radians.
        :returns: distance from station to the shower core in shower
                  front plane in m.

        """
        x = x - x0
        y = y - y0

        return sqrt(x ** 2 + y ** 2 -
                    (x * cos(phi) + y * sin(phi)) ** 2 * sin(theta) ** 2)


class NkgLdf(BaseLdf):

    """The Nishimura-Kamata-Greisen function"""

    # shower parameters
    # Age parameter and Moliere radius from Thoudam2012 sec 5.6.
    _Ne = 10 ** 4.8
    _s = 1.7
    _r0 = 30.

    def __init__(self, Ne=None, s=None):
        if Ne is not None:
            self._Ne = Ne
        if s is not None:
            self._s = s

        self._cache_c_s_value()

    def _cache_c_s_value(self):
        """Store the c_s value

        The c_s value does not change if s and r0 are fixed.

        """
        self._c_s = self._c(self._s)

    def calculate_ldf_value(self, r):
        """Calculate the LDF value for a given core distance

        :param r: core distance in m.
        :returns: particle density in m ** -2.

        """
        return self.get_ldf_value_for_size_and_shape(r, self._Ne, self._s)

    def get_ldf_value_for_size(self, r, Ne):
        """Calculate the LDF value for a given core distance and shower size

        :param r: core distance in m.
        :param Ne: number of electrons in the shower.
        :returns: particle density in m ** -2.

        """
        return self.get_ldf_value_for_size_and_shape(r, Ne, self._s)

    def get_ldf_value_for_size_and_shape(self, r, Ne, s):
        """Calculate the LDF value

        Given a core distance, shower size, and shower age.
        As given in Fokkema2012 eq 7.2.

        :param r: core distance in m.
        :param Ne: number of electrons in the shower.
        :param s: shower age parameter.
        :returns: particle density in m ** -2.

        """
        if s == self._s:
            c_s = self._c_s
        else:
            c_s = self._c(s)
        r0 = self._r0

        return Ne * c_s * (r / r0) ** (s - 2) * (1 + r / r0) ** (s - 4.5)

    def _c(self, s):
        """Part of the LDF

        As given in Fokkema2012 eq 7.3.

        :param s: shower age parameter.
        :returns: c(s)

        """
        r0 = self._r0
        return (gamma(4.5 - s) /
                (2 * pi * r0 ** 2 * gamma(s) * gamma(4.5 - 2 * s)))


class KascadeLdf(NkgLdf):

    """The KASCADE modified NKG function"""

    # shower parameters
    # Values from Fokkema2012 sec 7.1.
    _Ne = 10 ** 4.8
    _s = .94  # Shape parameter
    _r0 = 40.
    _alpha = 1.5
    _beta = 3.6

    def get_ldf_value_for_size_and_shape(self, r, Ne, s):
        """Calculate the LDF value

        Given a core distance, shower size, and shower age.
        As given in Fokkema2012 eq 7.4.

        :param r: core distance in m.
        :param Ne: number of electrons in the shower.
        :param s: shower shape parameter.
        :returns: particle density in m ** -2.

        """
        if s == self._s:
            c_s = self._c_s
        else:
            c_s = self._c(s)
        r0 = self._r0
        alpha = self._alpha
        beta = self._beta

        return Ne * c_s * (r / r0) ** (s - alpha) * (1 + r / r0) ** (s - beta)

    def _c(self, s):
        """Part of the LDF

        As given in Fokkema2012 eq 7.5.

        :param s: shower shape parameter.
        :returns: c(s)

        """
        r0 = self._r0
        beta = self._beta
        alpha = self._alpha
        return (gamma(beta - s) /
                (2 * pi * r0 ** 2 * gamma(s - alpha + 2) *
                 gamma(alpha + beta - 2 * s - 2)))
