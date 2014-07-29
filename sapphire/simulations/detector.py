"""Common HiSPARC station response simulations

These are some common simulations for HiSPARC detectors.

"""
from math import sqrt
import random

import numpy as np

from .base import BaseSimulation


class HiSPARCSimulation(BaseSimulation):

    def simulate_gps_uncertainty(self):
        """Simulate uncertainty from GPS receiver"""

        return np.random.normal(0, 4.5)

    def simulate_signal_transport_time(self, n=1):
        """ Simulate transport times of scintillation light to the PMT

        Generates random transit times within a given distribution and adds it
        to the times the particles passed the detector.

        :param n: number of times to simulate
        :returns: list of signal transport times

        """
        numbers = np.random.random(n)
        if n < 25:
            dt = [2.5507 + 2.39885 * x if x < 0.39377 else
                  1.56764 + 4.89536 * x for x in numbers]
        else:
            dt = np.where(numbers < 0.39377,
                          2.5507 + 2.39885 * numbers,
                          1.56764 + 4.89536 * numbers)
        return dt

    def simulate_detector_mips(self, particle_momenta):
        """Simulate the detector signal response for particles"""
        mips = sum(self.simulate_detector_mip(p) for p in particle_momenta)

        return mips

    def simulate_detector_mip(self, p):
        """Simulate the detector signal for one particle

        :param p: tuple with the x, y, z components of the particle momentum.

        """
        px, py, pz = p
        # determination of lepton angle of incidence
        costheta = abs(pz) / sqrt(px * px + py * py + pz * pz)

        # Simulation of convoluted distribution of electron and
        # muon energy losses with the scintillator response
        y = random.random()

        if y < 0.3394:
            mip = (0.48 + 0.8583 * sqrt(y)) / costheta
        elif y < 0.4344:
            mip = (0.73 + 0.7366 * y) / costheta
        elif y < 0.9041:
            mip = (1.7752 - 1.0336 * sqrt(0.9267 - y)) / costheta
        else:
            mip = (2.28 - 2.1316 * sqrt(1 - y)) / costheta

        return mip
