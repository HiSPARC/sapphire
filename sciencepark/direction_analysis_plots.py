import tables
from pylab import *

import master


def main():
    plot_sciencepark_cluster()

def plot_sciencepark_cluster():
    cluster = master.ScienceParkCluster(range(501, 507))

    figure()
    xl, yl = [], []
    for station in cluster.stations:
        for detector in station.detectors:
            x, y = detector.get_xy_coordinates()
            xl.append(x)
            yl.append(y)
            scatter(x, y, c='black', s=3)
    axis('equal')

    data = array([xl, yl]).T
    savetxt('plots/sciencepark-locations.txt', data)


if __name__ == '__main__':
    main()
