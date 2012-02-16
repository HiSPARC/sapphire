import tables
from itertools import combinations


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
                scatter(x, y, c='r', s=5, edgecolor='none')
            x, y, alpha = station.get_xyalpha_coordinates()
            scatter(x, y, c='orange', s=10, edgecolor='none')

    def plot_coincidence_twice(self, index=0, multiplicity=3):
        clf()
        subplot(121)
        self._do_plot_coincidence(index, multiplicity, use_detectors=False)
        subplot(122)
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
        title("Coincidence (%d-fold) #%d\n%s\n%s" % (multiplicity, index, method, type(self.solver).__name__))

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
        mylog = vectorize(lambda x: log10(x) if x > 0 else -999.)
        x = linspace(-400, 400, 100)
        y = linspace(-400, 400, 100)
        chi_squared = zeros((len(x), len(y)))
        for i in range(len(x)):
            for j in range(len(y)):
                chi_squared[j][i] = solver.calculate_chi_squared_for_xy(x[i], y[j])
        contourf(x, y, mylog(chi_squared), 50, cmap=cm.rainbow)
        colorbar()

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
        scatter(coincidence['x'], coincidence['y'], c='b', s=10)

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
        scatter(x, y, c='r', s=size)

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
        r = sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
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


class KascadeLdf(object):
    # shower parameters
    _Ne = 10 ** 4.8
    _s = .94
    _r0 = 40.
    _alpha = 1.5
    _beta = 3.6

    def __init__(self, *args, **kwargs):
        self.cache_c_s_value()

    def cache_c_s_value(self):
        self._c_s = self._c(self._s)

    def run(self):
        self.cache_c_s_value()

        super(KascadeLdfSimulation, self).run()

    def calculate_ldf_value(self, r):
        Ne = self._Ne
        s = self._s
        c_s = self._c_s
        r0 = self._r0
        alpha = self._alpha
        beta = self._beta

        return Ne * c_s * (r / r0) ** (s - alpha) * (1 + r / r0) ** (s - beta)

    def _c(self, s):
        r0 = self._r0
        beta = self._beta
        alpha = self._alpha
        return gamma(beta - s) / (2 * pi * r0 ** 2 * gamma(s - alpha + 2) * gamma(alpha + beta - 2 * s - 2))


def plot_minimal_and_overdetermined_solver(index, multiplicity, use_detectors=False, use_circles=False):
    ldf = KascadeLdf()

    if not use_circles:
        minimal_solver = CorePositionSolver(ldf)
        overdetermined_solver = OverdeterminedCorePositionSolver(ldf)
    else:
        minimal_solver = CorePositionCirclesSolver(ldf)
        overdetermined_solver = OverdeterminedCorePositionCirclesSolver(ldf)

    minimal_solver_reconstruction = CoreReconstruction(data, '/ldfsim', solver=minimal_solver)
    overdetermined_solver_reconstruction = CoreReconstruction(data, '/ldfsim', solver=overdetermined_solver)

    subplot(121)
    minimal_solver_reconstruction.plot_coincidence(index, multiplicity, use_detectors)
    subplot(122)
    overdetermined_solver_reconstruction.plot_coincidence(index, multiplicity, use_detectors)


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    c = CoreReconstruction(data, '/ldfsim')
