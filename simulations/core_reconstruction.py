import tables
from itertools import combinations

import numpy as np
import pylab as plt

from sapphire.simulations.ldf import KascadeLdf


DATAFILE = 'data.h5'


class CoreReconstruction(object):
    def __init__(self, data, simulation_path, solver=None):
        self.data = data
        self.simulation = data.getNode(simulation_path)
        self.cluster = self.simulation._v_attrs.cluster

        if solver is None:
            self.solver = CorePositionSolver(KascadeLdf())
        else:
            self.solver = solver

    def draw_cluster(self):
        for station in self.cluster.stations:
            for detector in station.detectors:
                x, y = detector.get_xy_coordinates()
                plt.scatter(x, y, c='r', s=5, edgecolor='none')
            x, y, alpha = station.get_xyalpha_coordinates()
            plt.scatter(x, y, c='orange', s=10, edgecolor='none')

    def plot_coincidence_twice(self, index=0, multiplicity=3):
        plt.clf()
        plt.subplot(121)
        self._do_plot_coincidence(index, multiplicity, use_detectors=False)
        plt.subplot(122)
        self._do_plot_coincidence(index, multiplicity, use_detectors=True)

    def plot_coincidence(self, index=0, multiplicity=3, use_detectors=False):
        self._do_plot_coincidence(index, multiplicity, use_detectors)

    def _do_plot_coincidence(self, index=0, multiplicity=3, use_detectors=False):
        coincidence = self.get_coincidence_with_multiplicity(index,
                                                             multiplicity)

        self._plot_chi_squared_on_map(coincidence, use_detectors)
        self._plot_coincidence_on_map(coincidence)

        for event in self.get_events_from_coincidence(coincidence):
            self._plot_event_on_map(event)

        if use_detectors:
            method = "individual detector signal"
        else:
            method = "station-averaged signal"
        plt.title("Coincidence (%d-fold) #%d\n%s\n%s" % (multiplicity, index, method, type(self.solver).__name__))

    def get_coincidence_with_multiplicity(self, index, multiplicity):
        coincidences = self.simulation.coincidences.read()
        sel = coincidences.compress(coincidences[:]['N'] >= multiplicity)
        return sel[index]

    def get_events_from_coincidence(self, coincidence):
        events = []
        id = coincidence['id']

        for index in self.simulation.c_index[id]:
            events.append(self.simulation.observables[index])

        return events

    def _plot_chi_squared_on_map(self, coincidence, use_detectors=False):
        solver = self.solver
        solver.reset_measurements()

        if use_detectors:
            self._add_detector_measurements_to_solver(solver, coincidence)
        else:
            self._add_station_measurements_to_solver(solver, coincidence)

        self._plot_chi_squared_contours(solver)

    def _plot_chi_squared_contours(self, solver):
        mylog = np.vectorize(lambda x: np.log10(x) if x > 0 else -999.)
        x = np.linspace(-400, 400, 100)
        y = np.linspace(-400, 400, 100)
        chi_squared = np.zeros((len(x), len(y)))
        for i in range(len(x)):
            for j in range(len(y)):
                chi_squared[j][i] = solver.calculate_chi_squared_for_xy(x[i], y[j])
        plt.contourf(x, y, mylog(chi_squared), 50, cmap=plt.cm.rainbow)
        plt.colorbar()

    def _add_station_measurements_to_solver(self, solver, coincidence):
        for event in self.get_events_from_coincidence(coincidence):
            if self._station_has_triggered(event):
                value = sum([event[u] for u in ['n1', 'n2', 'n3', 'n4']]) / 2.
                solver.add_measurement_at_xy(event['x'], event['y'], value)

    def _add_detector_measurements_to_solver(self, solver, coincidence):
        for event in self.get_events_from_coincidence(coincidence):
            station = self._get_station_from_event(event)
            if self._station_has_triggered(event):
                for detector, idx in zip(station.detectors, ['n1', 'n2', 'n3', 'n4']):
                    x, y = detector.get_xy_coordinates()
                    value = event[idx] / .5
                    solver.add_measurement_at_xy(x, y, value)

    def _plot_coincidence_on_map(self, coincidence):
        plt.scatter(coincidence['x'], coincidence['y'], c='b', s=10)

    def _plot_event_on_map(self, event):
        station_id = event['station_id']
        station = self.cluster.stations[station_id - 1]

        detectors = station.detectors
        num_particles = [event[u] for u in ['n1', 'n2', 'n3', 'n4']]

        for detector, num in zip(detectors, num_particles):
            self._plot_detector_on_map(detector, num)

    def _plot_detector_on_map(self, detector, num_particles):
        x, y = detector.get_xy_coordinates()
        size = num_particles * 10
        plt.scatter(x, y, c='r', s=size)

    def _station_has_triggered(self, event):
        if event['N'] >= 2:
            return True
        else:
            return False

    def _get_station_from_event(self, event):
        station_id = event['station_id']
        return self.cluster.stations[station_id - 1]


class CorePositionSolver(object):
    def __init__(self, ldf):
        self._measurements = []
        self._ldf = ldf

    def add_measurement_at_xy(self, x, y, value):
        self._measurements.append((x, y, value))

    def reset_measurements(self):
        self._measurements = []

    def calculate_chi_squared_for_xy(self, guess_x, guess_y):
        chi_squared = 0
        for expected, observed in self._get_expected_observed(guess_x, guess_y):
            chi_squared += (expected - observed) ** 2 / expected
        return chi_squared

    def _get_expected_observed(self, guess_x, guess_y):
        x0, y0, value0 = self._measurements[0]
        ldf_value0 = self._calculate_ldf_value_for_xy_xy(x0, y0, guess_x, guess_y)

        for x, y, value in self._measurements[1:]:
            ldf_value = self._calculate_ldf_value_for_xy_xy(x, y, guess_x, guess_y)
            expected = ldf_value / ldf_value0
            observed = value / value0
            yield expected, observed

    def _calculate_ldf_value_for_xy_xy(self, x0, y0, x1, y1):
        r = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
        return self._ldf.calculate_ldf_value(r)


class OverdeterminedCorePositionSolver(CorePositionSolver):
    def _get_expected_observed(self, guess_x, guess_y):
        for (x1, y1, value1), (x2, y2, value2) in combinations(self._measurements, 2):
            ldf_value1 = self._calculate_ldf_value_for_xy_xy(x1, y1, guess_x, guess_y)
            ldf_value2 = self._calculate_ldf_value_for_xy_xy(x2, y2, guess_x, guess_y)
            expected = ldf_value1 / ldf_value2
            observed = value1 / value2
            yield expected, observed


class CorePositionCirclesSolver(CorePositionSolver):
    def calculate_chi_squared_for_xy(self, guess_x, guess_y):
        chi_squared = 1
        for expected, observed in self._get_expected_observed(guess_x, guess_y):
            chi_squared *= (expected - observed) ** 2 / expected
        return chi_squared


class OverdeterminedCorePositionCirclesSolver(CorePositionCirclesSolver, OverdeterminedCorePositionSolver):
    pass


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    c = CoreReconstruction(data, '/ldfsim')
