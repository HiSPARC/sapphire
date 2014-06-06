"""
Lateral distribution functions that can be used for simulating particle
densities and for fitting to data.

"""
from scipy.special import gamma
from numpy import pi, sin, cos, sqrt


class BaseLdf(object):
    def calculate_core_distance_from_coordinates_and_direction(self,
                                                               x0, y0, x1, y1,
                                                               theta, phi):
        """Calculate core distance

        The core distance is the distance of the detector to the shower core,
        measured *on the shower front*.  For derivations, see logbook.

        """
        x = x1 - x0
        y = y1 - y0

        return sqrt(x ** 2 + y ** 2 -
                    (x * cos(phi) + y * sin(phi)) ** 2 * sin(theta) ** 2)


class KascadeLdf(BaseLdf):
    # shower parameters
    _Ne = 10 ** 4.8
    _s = .94
    _r0 = 40.
    _alpha = 1.5
    _beta = 3.6

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
