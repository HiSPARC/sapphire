"""Show Science Park detector locations on OpenStreetMap"""

import pylab as plt
import numpy as np
from matplotlib.path import Path
import matplotlib.patches as patches

import sapphire.api
import sapphire.clusters
import sapphire.simulations
from sapphire.simulations.ldf import KascadeLdf


DETECTOR_COLORS = ['black', 'r', 'g', 'b']
PATH_CODES = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO,
              Path.CLOSEPOLY]


def sciencepark_stations():
    try:
        network = sapphire.api.Network()
        stations = network.stations(subcluster=500)
        return [station['number'] for station in stations]
    except:
        return [501, 502, 503, 504, 505, 506, 508, 509]


def get_cluster(stations):
    cluster = sapphire.clusters.ScienceParkCluster(stations=stations)
    return cluster


def plot_detector_locations(cluster, background_path="backgrounds/ScienceParkMap_0.365.png"):
    plot_scintillators_in_cluster(cluster)
    draw_background_map(background_path)


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
            patch = patches.PathPatch(path, facecolor=DETECTOR_COLORS[i], lw=2)
            ax.add_patch(patch)
        station_x, station_y, station_a = station.get_xyalpha_coordinates()
        plt.scatter(station_x, station_y, marker='o', c='m', edgecolor='none',
                    s=7)
    plt.title('Science Park detector locations')
    plt.xlabel('Easting (meters)')
    plt.ylabel('Northing (meters)')


def draw_background_map(background_path):
    # Draw Science Park map on 1:1 scale (1 meter = 1 pixel)
    background = plt.imread(background_path)
    # determine pixel:meter ratio for different OSM zoom levels..
    bg_scale = 0.365
    bg_width = background.shape[1] * bg_scale
    bg_height = background.shape[0] * bg_scale
    plt.imshow(background, aspect='equal', alpha=0.5,
               extent=[-bg_width, bg_width, -bg_height, bg_height])


if __name__=="__main__":
    stations = sciencepark_stations()
    cluster = get_cluster(stations)
    plt.figure()
    plot_detector_locations(cluster)
    plt.show()
