import itertools

import numpy as np
import matplotlib.pyplot as plt

from sapphire.clusters import HiSPARCStations, ScienceParkCluster, SingleDiamondStation
from sapphire.analysis.direction_reconstruction import (DirectAlgorithm,
                                                        DirectAlgorithmCartesian3D,
                                                        DirectAlgorithmCartesian2D,
                                                        FitAlgorithm)


TIME_RESOLUTION = 2.5  # nanoseconds
C = .3  # lightspeed m/ns


def generate_discrete_times(station, detector_ids=[0, 2, 3]):
    """Generates possible arrival times for detectors

    The times are relative to the first detector, which is assumed to be
    at t = 0.

    """
    r = station_size(station, detector_ids)
    max_dt = ceil_in_base(r / C, TIME_RESOLUTION)
    times = np.arange(-max_dt, max_dt, TIME_RESOLUTION)
    time_combinations = itertools.product(times, repeat=len(detector_ids) - 1)
    return time_combinations


def station_size(station, detector_ids=[0, 2, 3]):
    """Get the largest distance between any two detectors in a station

    :param detectors: list of :class:`sapphire.clusters.Detector` objects

    """
    r = [station.calc_r_and_phi_for_detectors(d0, d1)[0]
         for d0, d1 in itertools.combinations(detector_ids, 2)]
    return max(r)


def ceil_in_base(value, base):
    return base * np.ceil(value / base)


if __name__ == '__main__':

    station_number = 502
    dirrec = DirectAlgorithmCartesian3D()

    try:
        station = HiSPARCStations([station_number]).get_station(station_number)
    except:
        station = ScienceParkCluster([station_number]).get_station(station_number)
    #station = SingleDiamondStation().stations[0]

    fig = plt.figure(figsize=(15, 10))
    sets = [plt.subplot2grid((2,3), (0,0), projection="polar"),
            plt.subplot2grid((2,3), (1,0), projection="polar"),
            plt.subplot2grid((2,3), (0,1), projection="polar"),
            plt.subplot2grid((2,3), (1,1), projection="polar")]
    combined = plt.subplot2grid((2,3), (1,2), projection="polar")
    layout = plt.subplot2grid((2,3), (0,2))

    # plt.setp(sets[0].get_xticklabels(), visible=False)
    # plt.setp(sets[2].get_xticklabels(), visible=False)
    # plt.setp(sets[2].get_yticklabels(), visible=False)
    # plt.setp(sets[3].get_yticklabels(), visible=False)
    # plt.setp(combined.get_yticklabels(), visible=False)

    detectors = [station.detectors[id].get_coordinates() for id in [0, 1, 2, 3]]
    x, y, z = zip(*detectors)

    layout.margins(0.25)
    layout.axis('equal')
    layout.scatter(x, y, s=15, marker='o', color='black')
    for id in [0, 1, 2, 3]:
        layout.annotate('%d' % id, (x[id], y[id]), xytext=(3, 3),
                        textcoords='offset points')
    layout.set_ylabel('northing (m)')
    layout.set_xlabel('easting (m)')

    colors = ['black', 'red', 'blue', 'green']

    for i, ids in enumerate(itertools.combinations([0, 1, 2, 3], 3)):
        times = generate_discrete_times(station, detector_ids=ids)
        detectors = [station.detectors[id].get_coordinates() for id in ids]
        x, y, z = zip(*detectors)

        theta, phi = itertools.izip(*(dirrec.reconstruct_common((0,) + t, x, y, z)
                                      for t in times))

        thetaa = np.degrees(np.array([t for t in theta if not np.isnan(t)]))
        phia = [p for p in phi if not np.isnan(p)]
        sets[i].scatter(phia, thetaa, s=1, marker='o', color='black')
        combined.scatter(phia, thetaa, s=1, marker='o', color=colors[i])

        sets[i].set_title(ids)
        sets[i].set_ylim(0, 90)
        sets[i].set_xlim(-185, 185)

    for i, ids in enumerate(itertools.combinations([0, 1, 2, 3], 3)):
        detectors = [station.detectors[id].get_coordinates() for id in ids]
        x, y, z = zip(*detectors)
        for t1 in (0, 10, 20, 30):
            times = ((t1, x) for x in np.arange(-60, 60, TIME_RESOLUTION))
            theta, phi = itertools.izip(*(dirrec.reconstruct_common((0,) + t, x, y, z)
                                          for t in times))
            thetaa = np.degrees(np.array([t for t in theta if not np.isnan(t)]))
            phia = [p for p in phi if not np.isnan(p)]
            sets[i].plot(phia, thetaa, color='red')

    combined.set_ylim(0, 90)
    combined.set_xlim(-185, 185)

    sets[0].set_ylabel('Zenith (degrees)')
    sets[3].set_xlabel('Azimuth (degrees)')
    fig.suptitle('Station: %d - Time resolution: %.1f ns' %
                      (station_number, TIME_RESOLUTION))
    plt.show()
