from detector_sim import *

for station in cluster:
    X, Y, beta = get_station_coordinates(station, 50, pi / 4, pi / 8)
    for detector in station['detectors']:
        x, y, orientation = detector
        c = get_detector_corners(X, Y, x, y, orientation, beta)
        cx, cy = zip(*c)
        fill(cx, cy)
