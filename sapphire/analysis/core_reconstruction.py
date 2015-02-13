""" Core reconstruction

    This module contains two classes that can be used to reconstruct
    HiSPARC events and coincidences. The classes know how to extract the
    relevant information from the station and event or cluster and
    coincidence. Various algorithms which do the reconstruction are also
    defined here. The algorithms require positions and particle densties to
    do the reconstruction.

    Each algorithm has a :meth:`~CenterMassAlgorithm.reconstruct_common`
    method which always requires particle denisties, x, and y positions
    and optionally z positions and previous reconstruction results. The
    data is then prepared for the algorithm and passed to the
    `reconstruct` method which returns the reconstructed x and y
    coordinates.

"""
from __future__ import division
import itertools

from numpy import isnan, nan, cos, sqrt, mean

from ..utils import pbar, ERR


class EventCoreReconstruction(object):

    """Reconstruct core for station events

    This class is aware of 'events' and 'stations'.  Initialize this class
    with a 'station' and you can reconstruct events using
    :meth:`reconstruct_event`.

    :param station: :class:`sapphire.clusters.Station` object.

    """

    def __init__(self, station):
        self.estimator = CenterMassAlgorithm
        self.station = station
        detectors = [d.get_coordinates() for d in self.station.detectors]
        self.area = [d.get_area() for d in self.station.detectors]
        self.x, self.y, self.z = zip(*detectors)

    def reconstruct_event(self, event, detector_ids=None):
        """Reconstruct a single event

        :param event: an event (e.g. from an events table), or any
            dictionary-like object containing the keys necessary for
            reconstructing the direction of a shower (e.g. number of mips).
        :param detector_ids: list of the detectors to use for
            reconstruction. The detector ids are 0-based, unlike the
            column names in the esd data.
        :returns: (x, y) core position in m.

        """
        p, x, y, z = ([], [], [], [])
        if detector_ids is None:
            detector_ids = range(4)
        for id in detector_ids:
            if event['n%d' % (id + 1)] not in ERR:
                p.append(event['n%d' % (id + 1)] / self.area[id])
                x.append(self.x[id])
                y.append(self.y[id])
                z.append(self.z[id])
        if len(p) >= 3:
            core_x, core_y = self.estimator.reconstruct_common(p, x, y, z)
        else:
            core_x, core_y = (nan, nan)
        return core_x, core_y

    def reconstruct_events(self, events, detector_ids=None, progress=True):
        """Reconstruct events

        :param events: the events table for the station from an ESD data
                       file.
        :param detector_ids: detectors which use for the reconstructions.
        :returns: (x, y) core positions in m.

        """
        cores = [self.reconstruct_event(event, detector_ids)
                 for event in pbar(events, show=progress)]
        core_x, core_y = zip(*cores)
        return core_x, core_y


class CoincidenceCoreReconstruction(object):

    """Reconstruct core for coincidences

    This class is aware of 'coincidences' and 'clusters'.  Initialize
    this class with a 'cluster' and you can reconstruct a coincidence
    using :meth:`reconstruct_coincidence`.

    :param cluster: :class:`sapphire.clusters.BaseCluster` object.

    """

    def __init__(self, cluster):
        self.estimator = CenterMassAlgorithm
        self.cluster = cluster

        # Store locations that do not change
        for station in cluster.stations:
            station.center_of_mass_coordinates = \
                station.calc_center_of_mass_coordinates()
            station.area = station.get_area()

    def reconstruct_coincidence(self, coincidence, station_numbers=None):
        """Reconstruct a single coincidence

        :param coincidence: a coincidence list consisting of
                            multiple (station_number, event) tuples
        :param station_numbers: list of station numbers, to only use
                                events from those stations.
        :returns: (x, y) core position in m.

        """
        p, x, y, z = ([], [], [], [])

        for station_number, event in coincidence:
            if station_numbers is not None:
                if station_number not in station_numbers:
                    continue
            try:
                sum_n = sum(event['n%d' % (i + 1)] for i in range(4)
                            if event['n%d' % (i + 1)] not in ERR)
            except ValueError:
                # All values -1 or -999
                continue
            station = self.cluster.get_station(station_number)
            p.append(sum_n / station.area)
            sx, sy, sz = station.center_of_mass_coordinates
            x.append(sx)
            y.append(sy)
            z.append(sz)

        if len(p) >= 3:
            core_x, core_y = self.estimator.reconstruct_common(p, x, y, z)
        else:
            core_x, core_y = (nan, nan)
        return core_x, core_y

    def reconstruct_coincidences(self, coincidences, station_numbers=None,
                                 progress=True):
        """Reconstruct all coincidences

        :param coincidences: a list of coincidences, each consisting of
                             multiple (station_number, event) tuples.
        :param station_numbers: list of station numbers, to only use
                                events from those stations.
        :returns: (x, y) core positions in m.

        """
        cores = [self.reconstruct_coincidence(coincidence, station_numbers)
                 for coincidence in pbar(coincidences, show=progress)]
        core_x, core_y = zip(*cores)
        return core_x, core_y


class CenterMassAlgorithm(object):

    """ Simple core estimator

    Estimates the core by center of mass of the measurements.

    """

    @classmethod
    def reconstruct_common(cls, p, x, y, z=None, initial={}):
        """Reconstruct core position

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param z: height of detectors is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions.

        """
        theta = initial.get('theta', nan)
        if not isnan(theta):
            p = [density * cos(theta) for density in p]

        return cls.reconstruct(p, x, y)

    @staticmethod
    def reconstruct(p, x, y):
        """Calculate center of mass

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.

        """
        core_x = sum(density * xi for density, xi in zip(p, x)) / sum(p)
        core_y = sum(density * yi for density, yi in zip(p, y)) / sum(p)
        return core_x, core_y


class AverageIntersectionAlgorithm(object):

    """ Core estimator

    To the densities in 3 stations correspond 2 possible cores. The line
    through these points is quit stable for the lateral distribution function.
    To each combination of 3 stations out of a set of at least 4
    stations hit corresponds a line. To each combinations of 2 lines out of
    the set of lines corresponds a point of intersection (if the 2 lines are
    not colinear). Taking the cloud of intersection points close to the core
    estimated by the center of mass, and averaging the positions in this cloud
    results in an estimation for the core.

    """

    @classmethod
    def reconstruct_common(cls, p, x, y, z=None, initial={}):
        """Reconstruct core

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param z: height of detectors is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions.

        """
        statindex = range(len(p))
        subsets = itertools.combinations(statindex, 3)

        m = 2.7  # average value in powerlaw  r ^(-m)  for density
        linelist0 = []
        linelist1 = []
        for zero, one, two in subsets:
            p0 = max(p[zero], .01)
            p1 = max(p[one], .01)
            p2 = max(p[two], .01)
            pp = (p0 / p1) ** (2. / m)
            qq = (p0 / p2) ** (2. / m)
            if pp == 1:
                pp = 1.000001
            if qq == 1:
                qq = 1.000001

            x0 = x[zero]
            x1 = x[one]
            x2 = x[two]
            y0 = y[zero]
            y1 = y[one]
            y2 = y[two]
            a = (x1 - pp * x0) / (1 - pp)
            b = (y1 - pp * y0) / (1 - pp)
            c = (x2 - qq * x0) / (1 - qq)
            d = (y2 - qq * y0) / (1 - qq)
            rsquare = pp * ((x1 - x0) ** 2 + (y1 - y0) ** 2) / ((1 - pp) ** 2)
            ssquare = qq * ((x2 - x0) ** 2 + (y2 - y0) ** 2) / ((1 - qq) ** 2)
            e = c - a
            f = d - b
            g = sqrt(e * e + f * f)
            k = 0.5 * (g * g + rsquare - ssquare) / g
            if abs(f) > 0:
                linelist0.append(-e / f)
                linelist1.append((a * e + b * f + g * k) / f)


        newx, newy = CenterMassAlgorithm.reconstruct_common(p, x, y, z,
                                                            initial)

        subsets = itertools.combinations(statindex, 2)

        xpointlist = []
        ypointlist = []
        for zero, one in subsets:
            a = linelist0[zero]
            b = linelist1[zero]
            c = linelist0[one]
            d = linelist1[one]
            aminc = max(a - c, 0.00001)
            xint = (d - b) / aminc
            yint = (a * d - b * c) / aminc

            if abs(xint - 2.5 * newx) < 200. and abs(yint - 2.5 * newy) < 200.:
                xpointlist.append(xint)
                ypointlist.append(yint)
        listmin = 1
        for distance in (100., 50., 25.):
            if len(xpointlist) > listmin:
                newx, newy, xpointlist, ypointlist = cls.select_newlist(
                    xpointlist, ypointlist, distance)
        if len(xpointlist) > listmin:
            newx = mean(xpointlist)
            newy = mean(ypointlist)

        return newx, newy

    @staticmethod
    def select_newlist(xpointlist, ypointlist, distance):
        """Select intersection points in square around the mean of old list."""

        newx = mean(xpointlist)
        newy = mean(ypointlist)
        newxlist = []
        newylist = []
        for xpoint, ypoint in zip(xpointlist, ypointlist):
            if abs(xpoint - newx) < distance and abs(ypoint - newy) < distance:
                newxlist.append(xpoint)
                newylist.append(ypoint)

        return newx, newy, newxlist, newylist
