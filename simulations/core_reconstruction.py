from __future__ import division

import tables
from itertools import combinations
import os
import sys

import numpy as np
import pylab as plt
from scipy import optimize
from scipy.misc import comb
from scipy.stats import scoreatpercentile
import progressbar as pb

from sapphire.simulations import ldf
from sapphire import storage

import utils

from pylab import *


DATAFILE = 'data.h5'

USE_TEX = False

# For matplotlib plots
if USE_TEX:
    rcParams['font.serif'] = 'Computer Modern'
    rcParams['font.sans-serif'] = 'Computer Modern'
    rcParams['font.family'] = 'sans-serif'
    rcParams['figure.figsize'] = [4 * x for x in (1, 2. / 3)]
    rcParams['figure.subplot.left'] = 0.175
    rcParams['figure.subplot.bottom'] = 0.175
    rcParams['font.size'] = 10
    rcParams['legend.fontsize'] = 'small'
    rcParams['text.usetex'] = True


class CoreReconstruction(object):
    def __init__(self, data, results_table=None, solver=None, N=None, overwrite=False):
        self.data = data

        if results_table:
            self.results_table = self.create_empty_output_table(results_table, overwrite)
        else:
            self.results_table = None

        if solver is None:
            self.solver = CorePositionSolver(ldf.KascadeLdf())
        else:
            self.solver = solver

        self.N = N

    def create_empty_output_table(self, table_path, overwrite=False):
        group, tablename = os.path.split(table_path)

        if table_path in self.data:
            if not overwrite:
                raise RuntimeError("Reconstruction table %s already exists" % table_path)
            else:
                self.data.removeNode(group, tablename)

        table = self._create_output_table(group, tablename)
        return table

    def _create_output_table(self, group, tablename):
        table = self.data.createTable(group, tablename,
                                      storage.ReconstructedEvent,
                                      createparents=True)
        return table

    def store_reconstructed_event(self, coincidence, event, reconstructed_core_x,
                                  reconstructed_core_y, reconstructed_shower_size):
        dst_row = self.results_table.row

        dst_row['id'] = event['id']
        dst_row['station_id'] = event['station_id']
        dst_row['r'] = coincidence['r']
        dst_row['phi'] = coincidence['phi']
        dst_row['alpha'] = event['alpha']
        dst_row['t1'] = event['t1']
        dst_row['t2'] = event['t2']
        dst_row['t3'] = event['t3']
        dst_row['t4'] = event['t4']
        dst_row['n1'] = event['n1']
        dst_row['n2'] = event['n2']
        dst_row['n3'] = event['n3']
        dst_row['n4'] = event['n4']
        dst_row['reference_theta'] = coincidence['shower_theta']
        dst_row['reference_phi'] = coincidence['shower_phi']
        dst_row['reference_core_pos'] = coincidence['x'], coincidence['y']
        dst_row['reconstructed_core_pos'] = reconstructed_core_x, reconstructed_core_y
        dst_row['reference_shower_size'] = coincidence['shower_size']
        dst_row['reconstructed_shower_size'] = reconstructed_shower_size
        dst_row['min_n134'] = min(event['n1'], event['n3'], event['n4'])
        dst_row.append()


    def reconstruct_core_positions(self, source):
        source = self.data.getNode(source)
        self.source = source

        self.cluster = source._v_attrs.cluster
        self._store_cluster_with_results()

        progressbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                              pb.ETA()],
                                     fd=sys.stderr)

        observables = source.observables
        coincidence_table = source.coincidences

        for event, coincidence in progressbar(zip(observables[:self.N], coincidence_table[:self.N])):
            assert event['id'] == coincidence['id']

            if coincidence['N'] >= 1:
                x, y, N = self.reconstruct_core_position(coincidence)
                self.store_reconstructed_event(coincidence, event, x, y, N)

        self.results_table.flush()

    def reconstruct_core_position(self, coincidence):
        solver = self.solver

        solver.reset_measurements()
        self._add_individual_detector_measurements_to_solver(solver, coincidence)

        x0, y0 = solver.get_center_of_mass_of_measurements()
        xopt, yopt = optimize.fmin(solver.calculate_chi_squared_for_xy, (x0, y0), disp=0)

        r, dens = solver.get_ldf_measurements_for_core_position((xopt, yopt))
        sigma = where(dens > 1, sqrt(dens), 1.)
        popt, pcov = optimize.curve_fit(solver.ldf_given_size, r, dens, p0=(1e5,), sigma=sigma)
        shower_size, = popt

        return xopt, yopt, shower_size

    def get_events_from_coincidence(self, coincidence):
        events = []
        id = coincidence['id']

        for index in self.source.c_index[id]:
            events.append(self.source.observables[index])

        return events

    def _add_station_averaged_measurements_to_solver(self, solver, coincidence):
        for event in self.get_events_from_coincidence(coincidence):
            if self._station_has_triggered(event):
                value = sum([event[u] for u in ['n1', 'n2', 'n3', 'n4']]) / 2.
                solver.add_measurement_at_xy(event['x'], event['y'], value)

    def _add_individual_detector_measurements_to_solver(self, solver, coincidence):
        for event in self.get_events_from_coincidence(coincidence):
            station = self._get_station_from_event(event)
            if self._station_has_triggered(event):
                for detector, idx in zip(station.detectors, ['n1', 'n2', 'n3', 'n4']):
                    x, y = detector.get_xy_coordinates()
                    value = event[idx] / .5
                    solver.add_measurement_at_xy(x, y, value)

    def _station_has_triggered(self, event):
        if event['N'] >= 2:
            return True
        else:
            return False

    def _get_station_from_event(self, event):
        station_id = event['station_id']
        return self.cluster.stations[station_id - 1]

    def _store_cluster_with_results(self):
        if not 'cluster' in self.results_table.attrs:
            self.results_table.attrs.cluster = self.cluster


class PlotCoreReconstruction(CoreReconstruction):
    def plot_reconstruct_core_position(self, source, coincidence_idx, multiplicity=1, min_detectors=3):
        source = self.data.getNode(source)
        self.source = source
        self.cluster = source._v_attrs.cluster
        coincidence = self.get_coincidence_with_multiplicity(coincidence_idx, multiplicity, min_detectors)

        figure()
        subplot(121)
        xopt, yopt, shower_size = self.reconstruct_core_position(coincidence)

        self._do_do_plot_coincidence(coincidence, use_detectors=True)
        x0, y0 = self.solver.get_center_of_mass_of_measurements()
        plt.scatter(x0, y0, color='green')
        plt.scatter(xopt, yopt, color='yellow')

        xlim(-50, 50)
        ylim(-50, 50)

        subplot(122)
        self._do_plot_ldf(coincidence, xopt, yopt, shower_size)

    def get_coincidence_with_multiplicity(self, index, multiplicity, num_detectors):
        coincidences = self.source.coincidences.read()
        sel = coincidences.compress(coincidences[:]['N'] >= multiplicity)
        coords = sel['id']

        events = self.source.observables.readCoordinates(coords)
        events_sel = events.compress(events[:]['N'] == num_detectors)

        idx = events_sel[index]['id']
        return coincidences[idx]

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

        self._do_do_plot_coincidence(coincidence, use_detectors)

    def _do_do_plot_coincidence(self, coincidence, use_detectors, index=0):
        self._plot_chi_squared_on_map(coincidence, use_detectors)
        self._plot_coincidence_on_map(coincidence)

        for event in self.get_events_from_coincidence(coincidence):
            self._plot_event_on_map(event)

        if use_detectors:
            method = "individual detector signal"
        else:
            method = "station-averaged signal"
        plt.title("Coincidence (%d-fold) #%d\n%s\n%s" % (coincidence['N'], index, method, type(self.solver).__name__))

    def _do_plot_ldf(self, coincidence, x0, y0, shower_size):
        x = logspace(-1, 2, 100)
        y = self.solver.ldf_given_size(x, shower_size)
        loglog(x, y, label='LDF')

        r, dens = self.solver.get_ldf_measurements_for_core_position((x0, y0))
        loglog(r, dens, 'o', label="signals")

        legend()
        xlabel("Core distance [m]")
        ylabel("Particle density [$m^{-2}$]")

    def _plot_chi_squared_on_map(self, coincidence, use_detectors=False):
        solver = self.solver
        solver.reset_measurements()

        if use_detectors:
            self._add_individual_detector_measurements_to_solver(solver, coincidence)
        else:
            self._add_station_measurements_to_solver(solver, coincidence)

        self._plot_chi_squared_contours(solver)

    def _plot_chi_squared_contours(self, solver):
        mylog = np.vectorize(lambda x: np.log10(x) if x > 0 else -999.)
        x = np.linspace(-100, 100, 100)
        y = np.linspace(-100, 100, 100)
        chi_squared = np.zeros((len(x), len(y)))
        for i in range(len(x)):
            for j in range(len(y)):
                chi_squared[j][i] = solver.calculate_chi_squared_for_xy((x[i], y[j]))
        plt.contourf(x, y, mylog(chi_squared), 100, cmap=plt.cm.rainbow)
        plt.colorbar()

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


class CorePositionSolver(object):
    _normalizing_idx = None

    def __init__(self, ldf):
        self._measurements = []
        self._ldf = ldf

    def add_measurement_at_xy(self, x, y, value):
        self._measurements.append((x, y, value))

    def reset_measurements(self):
        self._measurements = []
        self._normalizing_idx = None

    def calculate_chi_squared_for_xy(self, (guess_x, guess_y)):
        chi_squared = 0
        for expected, observed in self._get_expected_observed(guess_x, guess_y):
            if expected == np.Inf:
                return np.Inf
            else:
                chi_squared += (expected - observed) ** 2 / expected
        return chi_squared

    def get_ldf_measurements_for_core_position(self, (x0, y0)):
        r, ldf = [], []
        for x, y, value in self._measurements:
            core_distance = np.sqrt((x0 - x) ** 2 + (y0 - y) ** 2)
            if core_distance == 0.:
                core_distance = 1e-2
            r.append(core_distance)
            ldf.append(value)
        return array(r), array(ldf)

    def ldf_given_size(self, r, shower_size):
        return self._ldf.get_ldf_value_for_size(r, shower_size)

    def _get_expected_observed(self, guess_x, guess_y):
        normalizing_measurement, other_measurements = self._get_normalizing_and_other_measurements()

        x0, y0, value0 = normalizing_measurement
        ldf_value0 = self._calculate_ldf_value_for_xy_xy(x0, y0, guess_x, guess_y)

        for x, y, value in other_measurements:
            ldf_value = self._calculate_ldf_value_for_xy_xy(x, y, guess_x, guess_y)
            expected = ldf_value / ldf_value0
            observed = value / value0
            yield expected, observed

    def _get_normalizing_and_other_measurements(self):
        if self._normalizing_idx is not None:
            return self._split_out_measurement_at_idx(self._normalizing_idx)
        else:
            for idx, (x, y, value) in enumerate(self._measurements):
                if value > 0:
                    break

            self._normalizing_idx = idx

            if value <= 0:
                raise RuntimeError("BUG: serious problem with detector measurements!")
            else:
                return self._split_out_measurement_at_idx(idx)

    def _split_out_measurement_at_idx(self, idx):
        measurement = self._measurements[idx]
        other_measurements = self._measurements[:idx] + self._measurements[idx + 1:]
        return measurement, other_measurements

    def _calculate_ldf_value_for_xy_xy(self, x0, y0, x1, y1):
        r = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2)
        return self._ldf.calculate_ldf_value(r)

    def get_center_of_mass_of_measurements(self):
        measurements = np.array(self._measurements)
        x, y, value = measurements.T

        x0 = (value * x).sum() / value.sum()
        y0 = (value * y).sum() / value.sum()

        return x0, y0


class CorePositionSolverWithoutNullMeasurements(CorePositionSolver):
    def add_measurement_at_xy(self, x, y, value):
        if value > 0:
            self._measurements.append((x, y, value))


class OverdeterminedCorePositionSolver(CorePositionSolver):
    def _get_expected_observed(self, guess_x, guess_y):
        for (x1, y1, value1), (x2, y2, value2) in combinations(self._measurements, 2):
            ldf_value1 = self._calculate_ldf_value_for_xy_xy(x1, y1, guess_x, guess_y)
            ldf_value2 = self._calculate_ldf_value_for_xy_xy(x2, y2, guess_x, guess_y)
            expected = ldf_value1 / ldf_value2
            observed = value1 / value2
            yield expected, observed


class CorePositionCirclesSolver(CorePositionSolver):
    def calculate_chi_squared_for_xy(self, (guess_x, guess_y)):
        chi_squared = 1
        for expected, observed in self._get_expected_observed(guess_x, guess_y):
            chi_squared *= (expected - observed) ** 2 / expected
        return chi_squared

class CorePositionCirclesSolverWithoutNullMeasurements(CorePositionSolverWithoutNullMeasurements, CorePositionCirclesSolver):
    pass


class OverdeterminedCorePositionCirclesSolver(CorePositionCirclesSolver, OverdeterminedCorePositionSolver):
    pass


def do_reconstruction_plots(table):
    """Make plots based upon earlier reconstructions"""

    plot_N_reconstructions_vs_R(table)
    plot_core_pos_uncertainty_vs_R(table)
    plot_shower_size_hist(table)


Pnil = lambda x: exp(-0.5 * x)
Pp = lambda x: 1 - Pnil(x)
Ptrig = lambda x: comb(4, 2) * Pp(x) ** 2 * Pnil(x) ** 2 + \
                  comb(4, 3) * Pp(x) ** 3 * Pnil(x) + \
                  comb(4, 4) * Pp(x) ** 4

def plot_N_reconstructions_vs_R(table):
    figure()

    station = table.attrs.cluster.stations[0]

    x, y, alpha = station.get_xyalpha_coordinates()

    sim_path = table._v_pathname.replace('reconstructions', 'ldfsim')
    try:
        sim = data.getNode(sim_path)
    except tables.NoSuchNodeError:
        return

    # core distance for simulated events
    x2 = sim.coincidences.col('x')
    y2 = sim.coincidences.col('y')
    r = sqrt((x - x2) ** 2 + (y - y2) ** 2)

    # core distance for reconstructed events
    x2, y2 = table.col('reference_core_pos').T
    r2 = sqrt((x - x2) ** 2 + (y - y2) ** 2)

    bins = linspace(0, 50, 41)
    x, y = [], []
    for low, high in zip(bins[:-1], bins[1:]):
        sel = r.compress((low <= r) & (r < high))
        sel2 = r2.compress((low <= r2) & (r2 < high))

        if len(sel) > 0:
            x.append((low + high) / 2)
            y.append(len(sel2) / len(sel))
    x = array(x)
    y = array(y)

    plot(x, y, label="sim")

    kldf = ldf.KascadeLdf()
    dens = kldf.calculate_ldf_value(x)

    plot(x, Ptrig(dens), label="calc")
    legend()
    xlabel("Core distance [m]")
    ylabel("Reconstruction efficiency")
    utils.saveplot()

def plot_core_pos_uncertainty_vs_R(table):
    figure()

    x, y = table.col('reference_core_pos').T
    x2, y2 = table.col('reconstructed_core_pos').T
    d = sqrt((x - x2) ** 2 + (y - y2) ** 2)

    r = table.col('r')

    bins = linspace(0, 50, 41)
    x, d25, d50, d75 = [], [], [], []
    for low, high in zip(bins[:-1], bins[1:]):
        sel = d.compress((low <= r) & (r < high))

        if len(sel) > 0:
            x.append((low + high) / 2)
            d25.append(scoreatpercentile(sel, 25))
            d50.append(scoreatpercentile(sel, 50))
            d75.append(scoreatpercentile(sel, 75))

    fill_between(x, d25, d75, color='0.75')
    plot(x, d50, 'o-', color='black')

    xlabel("Core distance [m]")
    ylabel("Core position uncertainty [m]")
    utils.saveplot()

def plot_shower_size_hist(table):
    figure()

    reconstructed = table.col('reconstructed_shower_size')

    hist(log10(reconstructed), bins=200, histtype='step')
    axvline(log10(table[0]['reference_shower_size']))

    xlabel("log shower size")
    ylabel("count")
    utils.saveplot()


if __name__ == '__main__':
    np.seterr(divide='ignore')

    try:
        data
    except NameError:
        data = tables.openFile(DATAFILE, 'a')

    if '/reconstructions' not in data:
        c = CoreReconstruction(data, '/reconstructions/exact')
        c.reconstruct_core_positions('/ldfsim/exact')

        c = CoreReconstruction(data, '/reconstructions/gauss_10')
        c.reconstruct_core_positions('/ldfsim/gauss_10')

        c = CoreReconstruction(data, '/reconstructions/gauss_20')
        c.reconstruct_core_positions('/ldfsim/gauss_20')

        c = CoreReconstruction(data, '/reconstructions/poisson')
        c.reconstruct_core_positions('/ldfsim/poisson')

        c = CoreReconstruction(data, '/reconstructions/poisson_gauss_20')
        c.reconstruct_core_positions('/ldfsim/poisson_gauss_20')

        c = CoreReconstruction(data, '/reconstructions/poisson_gauss_20_nonull', solver=CorePositionSolverWithoutNullMeasurements(ldf.KascadeLdf()))
        c.reconstruct_core_positions('/ldfsim/poisson_gauss_20')

        c = CoreReconstruction(data, '/reconstructions/ground_gauss_20')
        c.reconstruct_core_positions('/groundsim/zenith_0/shower_0')

    utils.set_prefix("COR-")

    utils.set_suffix("-EXACT")
    do_reconstruction_plots(data.root.reconstructions.exact)

    utils.set_suffix("-GAUSS_10")
    do_reconstruction_plots(data.root.reconstructions.gauss_10)

    utils.set_suffix("-GAUSS_20")
    do_reconstruction_plots(data.root.reconstructions.gauss_20)

    utils.set_suffix("-POISSON")
    do_reconstruction_plots(data.root.reconstructions.poisson)

    utils.set_suffix("-POISSON-GAUSS_20")
    do_reconstruction_plots(data.root.reconstructions.poisson_gauss_20)

    utils.set_suffix("-POISSON-GAUSS_20_NONULL")
    do_reconstruction_plots(data.root.reconstructions.poisson_gauss_20_nonull)

    utils.set_suffix("-GROUND-GAUSS_20")
    do_reconstruction_plots(data.root.reconstructions.ground_gauss_20)
