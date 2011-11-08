import tables


DATAFILE = 'data.h5'


class CoreReconstruction(object):
    def __init__(self, data, simulation_path):
        self.data = data
        self.simulation = data.getNode(simulation_path)
        self.cluster = self.simulation._v_attrs.cluster

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
        self.do_plot_coincidence(index, multiplicity, use_detectors=False)
        subplot(122)
        self.do_plot_coincidence(index, multiplicity, use_detectors=True)

    def plot_coincidence(self, index=0, multiplicity=3, use_detectors=False):
        clf()
        self.do_plot_coincidence(index, multiplicity, use_detectors)

    def do_plot_coincidence(self, index=0, multiplicity=3, use_detectors=False):
        coincidence = self.get_coincidence_with_multiplicity(index,
                                                             multiplicity)

        self.plot_chi_squared_on_map(coincidence, use_detectors)
        self.plot_coincidence_on_map(coincidence)

        for event in self.get_events_from_coincidence(coincidence):
            self.plot_event_on_map(event)

        if use_detectors:
            method = "individual detector signal"
        else:
            method = "station-averaged signal"
        title("Coincidence (%d-fold) #%d\n%s" % (multiplicity, index, method))

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

    def plot_chi_squared_on_map(self, coincidence, use_detectors=False):
        solver = CorePositionSolver(KascadeLdf())

        if use_detectors:
            self.add_detector_measurements_to_solver(solver, coincidence)
        else:
            self.add_station_measurements_to_solver(solver, coincidence)

        self.plot_chi_squared_contours(solver)

    def plot_chi_squared_contours(self, solver):
        x = linspace(-400, 400, 100)
        y = linspace(-400, 400, 100)
        chi_squared = zeros((len(x), len(y)))
        for i in range(len(x)):
            for j in range(len(y)):
                chi_squared[j][i] = solver.calculate_chi_squared_for_xy(x[i], y[j])
        contourf(x, y, log10(chi_squared), 50, cmap=cm.rainbow)
        colorbar()

    def add_station_measurements_to_solver(self, solver, coincidence):
        for event in self.get_events_from_coincidence(coincidence):
            if self.station_has_triggered(event):
                value = sum([event[u] for u in ['n1', 'n2', 'n3', 'n4']]) / 2.
                solver.add_value_at_xy(event['x'], event['y'], value)

    def add_detector_measurements_to_solver(self, solver, coincidence):
        for event in self.get_events_from_coincidence(coincidence):
            station = self.get_station_from_event(event)
            if self.station_has_triggered(event):
                for detector, idx in zip(station.detectors, ['n1', 'n2', 'n3', 'n4']):
                    x, y = detector.get_xy_coordinates()
                    value = event[idx] / .5
                    solver.add_value_at_xy(x, y, value)

    def plot_coincidence_on_map(self, coincidence):
        scatter(coincidence['x'], coincidence['y'], c='b', s=10)

    def plot_event_on_map(self, event):
        station_id = event['station_id']
        station = self.cluster.stations[station_id - 1]

        detectors = station.detectors
        num_particles = [event[u] for u in ['n1', 'n2', 'n3', 'n4']]

        for detector, num in zip(detectors, num_particles):
            self.plot_detector_on_map(detector, num)

    def plot_detector_on_map(self, detector, num_particles):
        x, y = detector.get_xy_coordinates()
        size = num_particles * 10
        scatter(x, y, c='r', s=size)

    def station_has_triggered(self, event):
        if event['N'] >= 2:
            return True
        else:
            return False

    def get_station_from_event(self, event):
        station_id = event['station_id']
        return self.cluster.stations[station_id - 1]


class CorePositionSolver(object):
    def __init__(self, ldf):
        self.values = []
        self.ldf = ldf

    def add_value_at_xy(self, x, y, value):
        self.values.append((x, y, value))

    def calculate_chi_squared_for_xy(self, guess_x, guess_y):
        chi_squared = 0
        for expected, observed in self.get_expected_observed(guess_x, guess_y):
            chi_squared += (expected - observed) ** 2 / expected
        return chi_squared

    def get_expected_observed(self, guess_x, guess_y):
        x0, y0, value0 = self.values[0]
        ldf_value0 = self.calculate_ldf_value_for_xy_xy(x0, y0, guess_x, guess_y)

        for x, y, value in self.values[1:]:
            ldf_value = self.calculate_ldf_value_for_xy_xy(x, y, guess_x, guess_y)
            expected = ldf_value / ldf_value0
            observed = value / value0
            yield expected, observed

    def calculate_ldf_value_for_xy_xy(self, x0, y0, x1, y1):
        r = sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
        return self.ldf.calculate_ldf_value(r)


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


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'r')

    c = CoreReconstruction(data, '/ldfsim')
