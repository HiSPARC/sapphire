import tables
import numpy as np
import matplotlib.pyplot as plt

import utils


def main():
    global data
    data = tables.openFile('data-e15-S250.h5', 'r')
    utils.set_suffix('E_1PeV')

    scatterplot_core_distance_vs_time()


def scatterplot_core_distance_vs_time():
    figure()

    sim = data.root.showers.E_1PeV.zenith_0
    electrons = sim.electrons

    loglog(electrons[:]['core_distance'], electrons[:]['arrival_time'], ',')
    xlim(xmin=1e0)
    ylim(ymin=1e-3)

    xlabel("Core distance [m]")
    ylabel("Arrival time [ns]")

    utils.title("Shower front timing structure")
    utils.saveplot()


if __name__ == '__main__':
    main()
