""" This script is meant to illustrate the positioning of detectors
    in a station using the r, alpha, beta coordinates.

"""

from __future__ import division

import pylab as plt
import numpy as np
from math import pi, sqrt
from matplotlib.path import Path
import matplotlib.patches as patches

import sapphire.clusters
import sapphire.api

DETECTOR_COLORS = ['black', 'r', 'g', 'b']
PATH_CODES = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
              Path.CLOSEPOLY]


def get_cluster():
    cluster = sapphire.clusters.ClusterRAlphaBeta((0., 0.))
    # Diamond four detector station
    station_size = 10
    a = station_size / 2
    b = a * sqrt(3)
    c = sqrt(b ** 2 + station_size ** 2)
    p = np.arctan(2 / sqrt(3))
    q = pi / 2
    detectors = [(b, 0., 0.), (c, p, 0.),
                 (a, -q, q), (a, q, q)]
    cluster._add_station((-15, -15), detectors)
    # Triangle four detector station
    station_size = 10
    a = station_size / 2
    b = a * sqrt(3)
    q = pi / 2
    detectors = [(b, 0., 0.), (b / 3, 0., 0.),
                 (a, -q, q), (a, q, q)]
    cluster._add_station((15, -15), detectors)
    # Two detector station
    cluster._add_station((15, 15), [(5, -pi / 2, 0.),
                                    (5, pi / 2, 0.)])
    # Random - Swirl
#     detectors = [(2 + (i % 7), a, a)
#                  for i, a in enumerate(np.arange(0, 2 * pi, pi / 14))]
#     cluster._add_station((-15, 15), detectors)
    # HiSPARC Station from API (503, 504 or 505)
    station = sapphire.api.Station(504)
    detectors = station.detectors()
    detectors = [(d['radius'], np.radians(d['alpha']), np.radians(d['beta']))
                 for d in detectors]
    cluster._add_station((-15, 15), detectors)
    return cluster


def plot_scintillators_in_cluster(cluster):
    # Draw detector locations on a map
    ax = plt.gca()
    for station in cluster.stations:
        for i, detector in enumerate(station.detectors):
#             detector_x, detector_y = detector.get_xy_coordinates()
#             plt.scatter(detector_x, detector_y, marker='h',
#                         c=DETECTOR_COLORS[i], edgecolor='none', s=25)
            corners = detector.get_corners()
            corners.append(corners[0])  # Add first corner to complete outline
            path = Path(corners, PATH_CODES)
            j = i % len(DETECTOR_COLORS)
            patch = patches.PathPatch(path, facecolor=DETECTOR_COLORS[j], lw=2)
            ax.add_patch(patch)
        station_x, station_y, station_a = station.get_xyalpha_coordinates()
        plt.scatter(station_x, station_y, marker='o', c='m', edgecolor='none',
                    s=7)
    plt.title('Example r, alpha, beta detector locations')
    plt.xlabel('Easting (meters)')
    plt.ylabel('Northing (meters)')
    plt.axis('equal')
    plt.xlim(-40, 40)
    plt.ylim(-40, 40)

if __name__=="__main__":
    cluster = get_cluster()
    plt.figure()
    plot_scintillators_in_cluster(cluster)
    plt.show()
