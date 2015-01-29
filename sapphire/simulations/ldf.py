"""
Lateral distribution functions that can be used for simulating particle
densities and for fitting to data.

Example usage::

    import tables

    from sapphire.simulations.ldf import BaseLdfSimulation
    from sapphire.clusters import ScienceParkCluster

    data = tables.open_file('/tmp/test_ldf_simulation.h5', 'w')
    cluster = ScienceParkCluster()

    sim = BaseLdfSimulation(max_core_distance=400, cluster=cluster,
                            datafile=data, N=200)
    sim.run()

"""
from scipy.special import gamma
from numpy import pi, sin, cos, sqrt, random

from .detector import HiSPARCSimulation, ErrorlessSimulation
from ..utils import pbar


class BaseLdfSimulation(HiSPARCSimulation):

    def __init__(self, max_core_distance, *args, **kwargs):
        super(BaseLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = NkgLdf()
        self.max_core_distance = max_core_distance

        # Since the cluster is not rotated detector positions can be stored.
        for station in self.cluster.stations:
            for detector in station.detectors:
                detector.xy_coordinates = detector.get_xy_coordinates()

    def generate_shower_parameters(self):
        """Generate shower parameters, i.e. core position.

        For the simple LDF only the core position is relevant. It
        assumes the shower to come from the Zenith.

        :returns: dictionary with shower parameters: core_pos
                  (x, y-tuple).

        """
        r = self.max_core_distance
        giga = int(1e9)

        for i in pbar(range(self.N)):
            shower_parameters = {'ext_timestamp': (giga + i) * giga,
                                 'azimuth': 0.,
                                 'zenith': 0.,
                                 'core_pos': self.generate_core_position(r),
                                 'size': None,
                                 'energy': None}

            yield shower_parameters

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Get the mips in a detector from the LDF.

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
        x, y = detector.xy_coordinates
        core_x, core_y = shower_parameters['core_pos']
        zenith = shower_parameters['zenith']
        azimuth = shower_parameters['azimuth']

        r = self.ldf.calculate_core_distance_from_coordinates_and_direction(
                x, y, core_x, core_y, zenith, azimuth)
        p = self.ldf.calculate_ldf_value(r)
        num_particles = random.poisson(p *  detector.get_area())

        return num_particles


class KascadeLdfSimulation(BaseLdfSimulation):

    def __init__(self, *args, **kwargs):
        super(KascadeLdfSimulation, self).__init__(*args, **kwargs)

        self.ldf = KascadeLdf()


class BaseLdf(object):

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

    # shower parameters
    # Age parameter (shape) and Moliere radius from Thoudam2012 sec 5.6.
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
        self._c_s = self._c(self._s)

    def calculate_ldf_value(self, r):
        return self.get_ldf_value_for_size_and_shape(r, self._Ne, self._s)

    def get_ldf_value_for_size(self, r, Ne):
        return self.get_ldf_value_for_size_and_shape(r, Ne, self._s)

    def get_ldf_value_for_size_and_shape(self, r, Ne, s):
        c_s = self._c_s
        r0 = self._r0

        return Ne * c_s * (r / r0) ** (s - 2) * (1 + r / r0) ** (s - 4.5)

    def _c(self, s):
        r0 = self._r0
        return (gamma(4.5 - s) /
                (2 * pi * r0 ** 2 * gamma(s) * gamma(4.5 - 2 * s)))


class KascadeLdf(NkgLdf):

    # shower parameters
    # Values from Fokkema2012 sec 7.1.
    _Ne = 10 ** 4.8
    _s = .94
    _r0 = 40.
    _alpha = 1.5
    _beta = 3.6

    def get_ldf_value_for_size_and_shape(self, r, Ne, s):
        c_s = self._c_s
        r0 = self._r0
        alpha = self._alpha
        beta = self._beta

        return Ne * c_s * (r / r0) ** (s - alpha) * (1 + r / r0) ** (s - beta)

    def _c(self, s):
        r0 = self._r0
        beta = self._beta
        alpha = self._alpha
        return (gamma(beta - s) /
                (2 * pi * r0**2 * gamma(s - alpha + 2) *
                 gamma(alpha + beta - 2 * s - 2)))
