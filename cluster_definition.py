from __future__ import division

from math import sqrt, pi


DETECTOR_SIZE = (.5, 1.)
STATION_SIZE = 10
CLUSTER_SIZE = 250

a = STATION_SIZE / 2
b = a / 3 * sqrt(3)
detectors = [(0., 2 * b, 'UD'), (0., 0., 'UD'),
             (-a, -b, 'LR'), (a, -b, 'LR')]

A = CLUSTER_SIZE / 2
B = A / 3 * sqrt(3)
cluster = [{'position': (0, 0), 'angle': 0, 'detectors': detectors},
           {'position': (0, 2 * B), 'angle': 0, 'detectors': detectors},
           {'position': (-A, -B), 'angle': 2 * pi / 3,
            'detectors': detectors},
           {'position': (A, -B), 'angle': -2 * pi / 3,
            'detectors': detectors},
          ]
