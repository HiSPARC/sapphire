import numpy as np
import pylab as plt

from sapphire import clusters
from sapphire.simulations.ldf import KascadeLdf


class ToyMonteCarlo(object):
    def main(self):
        global weights, densities, weighted_densities
        plt.figure()

        cluster = clusters.SingleStation()
        self.station = cluster.stations[0]

        R = np.linspace(0, 100, 100)
        densities = []
        weights = []
        for E in np.linspace(1e13, 1e17, 10000):
            relative_flux = E ** -2.7
            Ne = 10 ** (np.log10(E) - 15 + 4.8)
            self.ldf = KascadeLdf(Ne)
            min_dens = self.calculate_minimum_density_for_station_at_R(R)

            weights.append(relative_flux)
            densities.append(min_dens)
        weights = np.array(weights)
        densities = np.array(densities).T

        weighted_densities = (np.sum(weights * densities, axis=1) /
                              np.sum(weights))
        plt.plot(R, weighted_densities)
        plt.yscale('log')
        plt.ylabel("Min. density [m^{-2}]")
        plt.xlabel("Core distance [m]")
        plt.axvline(5.77)
        plt.show()

    def calculate_minimum_density_for_station_at_R(self, R):
        densities = self.calculate_densities_for_station_at_R(R)
        return np.min(densities, axis=0)

    def calculate_densities_for_station_at_R(self, R):
        densities = []
        for detector in self.station.detectors:
            densities.append(self.calculate_densities_for_detector_at_R(
                                detector, R))
        return np.array(densities)

    def calculate_densities_for_detector_at_R(self, detector, R):
        x = 0
        y = R
        x0, y0 = detector.get_xy_coordinates()

        r = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
        return self.ldf.calculate_ldf_value(r)


if __name__ == '__main__':
    sim = ToyMonteCarlo()
    sim.main()
