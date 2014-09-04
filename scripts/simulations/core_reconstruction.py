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

from sapphire.simulations import ldf
from sapphire import storage
from sapphire.analysis.core_reconstruction import *

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


def do_reconstruction_plots(table):
    """Make plots based upon earlier reconstructions"""

    plot_N_reconstructions_vs_R(table)
    plot_core_pos_uncertainty_vs_R(table)
    plot_shower_size_hist(table)
    plot_scatter_reconstructed_core(table, N=1000)


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
        sim = data.get_node(sim_path)
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
    reference_shower_size = table[0]['reference_shower_size']
    if reference_shower_size == 0.:
        reference_shower_size = 10 ** 4.8
    axvline(log10(reference_shower_size))

    xlabel("log shower size")
    ylabel("count")
    utils.saveplot()


def plot_scatter_reconstructed_core(table, N=None):
    # Make sure to get a *copy*
    figsize = list(rcParams['figure.figsize'])
    figsize[0] = figsize[1] * 2

    figure(figsize=figsize)

    station = table.attrs.cluster.stations[0]
    subplot(121)
    x, y = table.col('reference_core_pos')[:N].T
    #scatter(x, y, c='b', s=1, edgecolor='none', zorder=1)
    plot(x, y, ',', c='b', markeredgecolor='b', zorder=1)
    for detector in station.detectors:
        x, y = detector.get_xy_coordinates()
        plt.scatter(x, y, c='r', s=20, edgecolor='none', zorder=2)
    xlabel("Distance [m]")
    ylabel("Distance [m]")
    xlim(-60, 60)
    ylim(-60, 60)
    title("simulated")

    subplot(122)
    x, y = table.col('reconstructed_core_pos')[:N].T
    #scatter(x, y, c='b', s=1, edgecolor='none', zorder=1)
    plot(x, y, ',', c='b', markeredgecolor='b', zorder=1)
    for detector in station.detectors:
        x, y = detector.get_xy_coordinates()
        plt.scatter(x, y, c='r', s=20, edgecolor='none', zorder=2)
    xlabel("Distance [m]")
    ylabel("Distance [m]")
    xlim(-60, 60)
    ylim(-60, 60)
    title("reconstructed")

    utils.saveplot()


if __name__ == '__main__':
    np.seterr(divide='ignore')

    try:
        data
    except NameError:
        data = tables.open_file(DATAFILE, 'a')

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

        #c = CoreReconstruction(data, '/reconstructions/ground_gauss_20')
        #c.reconstruct_core_positions('/groundsim/zenith_0/shower_0')

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

    #utils.set_suffix("-GROUND-GAUSS_20")
    #do_reconstruction_plots(data.root.reconstructions.ground_gauss_20)
