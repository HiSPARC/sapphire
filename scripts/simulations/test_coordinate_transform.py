import tables
from pylab import *


def plot_station_and_shower_transforms(event_id):
    # plot old coordinates
    coincidence = sim.coincidences[event_id]
    r, phi, alpha = coincidence['r'], coincidence['phi'], \
                    coincidence['alpha']
    cluster.set_rphialpha_coordinates(r, phi, alpha)
    plot_cluster(0.2)
    scatter(0, 0, c='r', alpha=.2)

    # plot new coordinates
    coincidence = test_output.coincidences[event_id]
    x, y = coincidence['x'], coincidence['y']
    cluster.set_rphialpha_coordinates(0, 0, 0)
    plot_cluster(1.)
    scatter(x, y, c='r', alpha=1.)
    scatter(0, 0, c='white', alpha=1.)
    scatter(0, 0, c='r', alpha=.2)

    # plot coordinates stored in 'observables' table
    for event in test_output.observables.read_where('id == %d' % event_id):
        scatter(event['x'], event['y'], c='lightgreen')

    xlabel("[m]")
    ylabel("[m]")


def plot_cluster(alpha):
    for station in cluster.stations:
        c = 'yellow'
        for detector in station.detectors:
            x, y = detector.get_xy_coordinates()
            scatter(x, y, c=c, alpha=alpha)
            if c == 'yellow':
                c = 'blue'


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.open_file('transform_test.h5')
        sim = data.root.simulations.E_100TeV.zenith_0
        test_output = data.root.test_output.E_100TeV.zenith_0
        cluster = test_output._v_attrs.cluster
