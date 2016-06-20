"""
Lateral distribution functions that can be used for simulating particle
densities and for fitting to data.

Example usage::

    >>> import tables

    >>> from sapphire import NkgLdfSimulation, ScienceParkCluster

    >>> data = tables.open_file('/tmp/test_ldf_simulation.h5', 'w')
    >>> cluster = ScienceParkCluster()

    >>> sim = NkgLdfSimulation(max_core_distance=400, min_energy=1e15,
    ...                        max_energy=1e21, cluster=cluster,
    ...                        datafile=data, N=200)
    >>> sim.run()

"""
import warnings

from scipy.special import gamma
from numpy import pi, sin, cos, sqrt, random, arctan2, log10

from .detector import HiSPARCSimulation, ErrorlessSimulation
from ..utils import pbar, vector_length, memoize


class BaseLdfSimulation(HiSPARCSimulation):

    def __init__(self, max_core_distance, min_energy, max_energy, *args,
                 **kwargs):
        """Simulation initialization

        :param max_core_distance: maximum distance of shower core to
                                  center of cluster (in meters).
        :param min_energy,max_energy: Minimum and maximum energy of the
                                      shower (in eV).

        """
        super(BaseLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = BaseLdf()
        self.max_core_distance = max_core_distance
        self.min_energy = min_energy
        self.max_energy = max_energy

        # The cluster is not moved, so detector positions can be stored.
        for station in self.cluster.stations:
            for detector in station.detectors:
                detector.xy_coordinates = detector.get_xy_coordinates()

    def generate_shower_parameters(self):
        """Generate shower parameters, i.e. core position

        For the simple LDF only the core position is relevant. It
        assumes the shower to come from the Zenith.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple).

        """
        r = self.max_core_distance
        giga = int(1e9)

        for i in pbar(range(self.N), show=self.progress):
            energy = self.generate_energy(self.min_energy, self.max_energy)
            size = 10 ** (log10(energy) - 15 + 4.8)
            shower_parameters = {'ext_timestamp': (giga + i) * giga,
                                 'azimuth': self.generate_azimuth(),
                                 'zenith': 0.,
                                 'core_pos': self.generate_core_position(r),
                                 'size': size,
                                 'energy': energy}

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
        size = shower_parameters['size']

        r = self.ldf.calculate_core_distance(x, y, core_x, core_y, zenith,
                                             azimuth)

        p_shower = self.ldf.calculate_ldf_value(r, Ne=size)
        p_ground = p_shower * cos(zenith)
        num_particles = self.simulate_particles_for_density(
            p_ground * detector.get_area())

        return num_particles

    @staticmethod
    def simulate_particles_for_density(p):
        """Get number of particles in detector given a particle density

        :param p: particle density in number per detector area.
        :return: random number from Poisson distribution.

        """
        return random.poisson(p)


class BaseLdfSimulationWithoutErrors(ErrorlessSimulation, BaseLdfSimulation):

    """This simulation does not simulate errors/uncertainties

    This should result in perfect particle counting for the detectors.

    """

    @staticmethod
    def simulate_particles_for_density(p):
        """Exact number"""

        return p


class NkgLdfSimulation(BaseLdfSimulation):

    """Same as the BaseLdfSimulation but uses the NkgLdf as LDF"""

    def __init__(self, *args, **kwargs):
        super(NkgLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = NkgLdf()


class NkgLdfSimulationWithoutErrors(NkgLdfSimulation,
                                    BaseLdfSimulationWithoutErrors):

    """Same as the NkgLdfSimulation but without error simulation"""

    pass


class KascadeLdfSimulation(BaseLdfSimulation):

    """Same as the BaseLdfSimulation but uses the KascadeLdf as LDF"""

    def __init__(self, *args, **kwargs):
        super(KascadeLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = KascadeLdf()


class KascadeLdfSimulationWithoutErrors(KascadeLdfSimulation,
                                        BaseLdfSimulationWithoutErrors):

    """Same as the KascadeLdfSimulation but without error simulation"""

    pass


class EllipseLdfSimulation(BaseLdfSimulation):

    """Same as BaseLdfSimulation but uses the EllipseLdF as LDF"""

    def __init__(self, *args, **kwargs):
        super(EllipseLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = EllipseLdf()

    def generate_shower_parameters(self):
        """Generate shower parameters, i.e. core position

        For the elliptic LDF both the core position and the zenith angle
        are relevant.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple).

        """
        r = self.max_core_distance
        giga = int(1e9)

        for i in pbar(range(self.N), show=self.progress):
            energy = self.generate_energy(self.min_energy, self.max_energy)
            size = 10 ** (log10(energy) - 15 + 4.8)
            shower_parameters = {'ext_timestamp': (giga + i) * giga,
                                 'azimuth': self.generate_azimuth(),
                                 'zenith': self.generate_zenith(),
                                 'core_pos': self.generate_core_position(r),
                                 'size': size,
                                 'energy': energy}

            yield shower_parameters

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
        size = shower_parameters['size']

        r, phi = self.ldf.calculate_core_distance_and_angle(x, y, core_x,
                                                            core_y)

        p_ground = self.ldf.calculate_ldf_value(r, phi, size, zenith, azimuth)
        num_particles = self.simulate_particles_for_density(
            p_ground * detector.get_area())

        return num_particles


class BaseLdf(object):

    def calculate_ldf_value(self, r, Ne=None, s=None):
        return 0.

    def calculate_core_distance(self, x, y, x0, y0, theta, phi):
        """Calculate core distance

        The core distance is the distance of the detector to the shower core,
        measured *on the shower front*.  For derivations, see logbook.

        :param x,y: detector position in m.
        :param x0,y0: shower core position in m.
        :param theta,phi: shower axis direction in radians.
        :return: distance from detector to the shower core in shower
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
        """NKG LDF setup

        :param Ne: Shower size (number of electrons).
        :param s: Shower age parameter.

        """
        if Ne is not None:
            self._Ne = Ne
        if s is not None:
            self._s = s

    def calculate_ldf_value(self, r, Ne=None, s=None):
        """Calculate the LDF value

        :param r: core distance in m.
        :param Ne: number of electrons in the shower.
        :param s: shower age parameter.
        :return: particle density in m^-2.

        """
        if Ne is None:
            Ne = self._Ne
        if s is None:
            s = self._s
        return self.ldf_value(r, Ne, s)

    def ldf_value(self, r, Ne, s):
        """Calculate the LDF value

        Given a core distance, shower size, and shower age.
        As given in Fokkema2012 eq 7.2.

        :param r: core distance in m.
        :param Ne: number of electrons in the shower.
        :param s: shower age parameter.
        :return: particle density in m^-2.

        """
        r0 = self._r0
        c_s = self._c(s, r0)
        return Ne * c_s * (r / r0) ** (s - 2) * (1 + r / r0) ** (s - 4.5)

    @memoize
    def _c(self, s, r0):
        """Part of the LDF

        As given in Fokkema2012 eq 7.3.

        :param s: shower age parameter.
        :param r0: Moliere radius.
        :return: c(s)

        """
        return (gamma(4.5 - s) /
                (2 * pi * r0 ** 2 * gamma(s) * gamma(4.5 - 2 * s)))


class KascadeLdf(NkgLdf):

    """The KASCADE modified NKG function"""

    # shower parameters
    # Values from Fokkema2012 sec 7.1.
    _Ne = 10 ** 4.8
    _s = 0.94  # Shape parameter
    _r0 = 40.
    _alpha = 1.5
    _beta = 3.6

    def ldf_value(self, r, Ne, s):
        """Calculate the LDF value

        Given a core distance, shower size, and shower age.
        As given in Fokkema2012 eq 7.4.

        :param r: core distance in m.
        :param Ne: number of electrons in the shower.
        :param s: shower shape parameter.
        :return: particle density in m^-2.

        """
        r0 = self._r0
        alpha = self._alpha
        beta = self._beta
        c_s = self._c(s, r0, alpha, beta)

        return Ne * c_s * (r / r0) ** (s - alpha) * (1 + r / r0) ** (s - beta)

    @memoize
    def _c(self, s, r0, alpha, beta):
        """Part of the LDF

        As given in Fokkema2012 eq 7.5.

        :param s: shower shape parameter.
        :param r0: Moliere radius.
        :param alpha: shape.
        :param beta: shape.
        :return: c(s)

        """
        return (gamma(beta - s) /
                (2 * pi * r0 ** 2 * gamma(s - alpha + 2) *
                 gamma(alpha + beta - 2 * s - 2)))


class EllipseLdf(KascadeLdf):

    """NKG function modified for electrons and muons and azimuthal asymmetry"""

    # shower parameters
    # Values from Montanus, paper to follow.
    _Ne = 10 ** 6.0   # Ne is number of electrons and muons
    _s1 = -0.592  # Shape parameter
    _s2 = -3.157  # Shape parameter
    _r0 = 30.
    _zenith = 0.
    _azimuth = 0.

    def __init__(self, Ne=None, zenith=None, azimuth=None, s1=None, s2=None):
        if Ne is not None:
            self._Ne = Ne
        if zenith is not None:
            self._zenith = zenith
        if azimuth is not None:
            self._azimuth = azimuth
        if s1 is not None:
            self._s1 = s1
        if s2 is not None:
            self._s2 = s2

    def calculate_ldf_value(self, r, phi, Ne=None, zenith=None, azimuth=None):
        """Calculate LDF value for given distance and polar angle

        :param r: core distance in m.
        :param phi: polar angle in rad.
        :param Ne: number of electrons and muons in the shower.
        :return: particle density in m^-2.

        """
        if Ne is None:
            Ne = self._Ne
        if zenith is None:
            zenith = self._zenith
        if azimuth is None:
            azimuth = self._azimuth
        return self.ldf_value(r, phi, Ne, zenith, azimuth, self._s1, self._s2)

    def ldf_value(self, r, phi, Ne, zenith, azimuth, s1, s2):
        """Calculate the LDF value for a given parameters

        Method as given by Montanus, paper to follow.

        :param r: core distance in m.
        :param phi: polar angle in rad.
        :param Ne: number of electrons and muons in the shower.
        :param zenith: zenith angle in rad.
        :param azimuth: azimuth angle in rad.
        :param s1: shower shape parameter.
        :param s2: shower shape parameter.
        :return: particle density in m^-2.

        """
        s1 = s1 + 0.229 * zenith
        s2 = s2 + 0.222 * zenith
        r0 = self._r0
        c_s = self._c(s1, s2, r0)

        # zenith = self._zenith
        # azimuth = self._azimuth
        # relative polar angle
        relpolar = cos(phi - azimuth)

        # elliptic iso-density contours
        ell = sqrt(1. - sin(zenith) ** 2 * relpolar ** 2)
        # shift of the center of elliptic iso-density contours
        shift = -0.058 * sin(2. * zenith) * relpolar

        k = r * (shift + ell)
        term1 = 1. * k / r0
        term2 = 1. + term1
        # empirical modification
        greis = 11.4 * cos(zenith) ** 2
        # c(s) * normcorr = normalization
        normcorr = greis / (2. + s1 - greis * (3 + s1 + s2))
        # Greisen modification of NKG LDF
        term3 = 1. + term1 / greis

        with warnings.catch_warnings(record=True):
            density = (Ne * cos(zenith) * c_s * term1 ** s1 * term2 ** s2 *
                       normcorr * term3)
        return density

    @memoize
    def _c(self, s1, s2, r0):
        """Normalization of the LDF

        As given in Montanus, paper to follow.

        :param s1: shower shape parameter.
        :param s2: shower shape parameter.
        :return: c(s1,s2).

        """
        return (gamma(-s2) /
                (2 * pi * r0 ** 2 * gamma(s1 + 2) * gamma(-s1 - s2 - 3)))

    def calculate_core_distance_and_angle(self, x, y, x0, y0):
        """Calculate core distance

        The core distance is the distance of the detector to the shower core,
        measured *in the horizontal observation plane*.

        :param x,y: detector position in m.
        :param x0,y0: shower core position in m.
        :return: distance in m, and polar angle in radians from shower core
                 to detector in horizontal observation plane.

        """
        dx = x - x0
        dy = y - y0

        return vector_length(dx, dy), arctan2(dy, dx)
