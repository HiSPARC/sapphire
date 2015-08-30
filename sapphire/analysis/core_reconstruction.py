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
    data is then prepared for the algorithm and passed to
    the :meth:`~CenterMassAlgorithm.reconstruct` method which returns the
    reconstructed x and y coordinates.

"""
from __future__ import division
import itertools

from numpy import isnan, nan, cos, sqrt, log10, mean, array

from .event_utils import station_density, detector_density
from ..utils import pbar
from ..simulations import ldf


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
        :return: (x, y) core position in m.

        """
        p, x, y, z = ([], [], [], [])
        if detector_ids is None:
            detector_ids = range(4)
        for id in detector_ids:
            p_detector = detector_density(event, id, self.station)
            if not isnan(p_detector):
                p.append(p_detector)
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
        :return: (x, y) core positions in m.

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
        :return: (x, y) core position in m.

        """
        p, x, y, z = ([], [], [], [])

        for station_number, event in coincidence:
            if station_numbers is not None:
                if station_number not in station_numbers:
                    continue
            station = self.cluster.get_station(station_number)
            p_station = station_density(event, range(4), station)
            if not isnan(p_station):
                sx, sy, sz = station.center_of_mass_coordinates
                p.append(p_station)
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
        :return: (x, y) core positions in m.

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
        #theta = initial.get('theta', nan)
        #if not isnan(theta):
            #p = [density * cos(theta) for density in p]

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
    not collinear). Taking the cloud of intersection points close to the core
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
        if len(p) < 4 or len(x) < 4 or len(y) < 4:
            raise Exception('This algorithm requires at least 4 detections.')
        phit = []
        xhit = []
        yhit = []
        for i in range(len(p)):
            if p[i] > .00001:
                phit.append(p[i])
                xhit.append(x[i])
                yhit.append(y[i])

        statindex = range(len(phit))
        subsets = itertools.combinations(statindex, 3)      # select triples of stations
        m = 2.3  # optimized value in powerlaw  r ^(-m)  for density

        linelist0 = []
        linelist1 = []
        for zero, one, two in subsets:
            pp = (phit[zero] / phit[one]) ** (2. / m)
            qq = (phit[zero] / phit[two]) ** (2. / m)
            if pp == 1:
                pp = 1.000001       # to avoid singularity
            if qq == 1:
                qq = 1.000001       # to avoid singularity

            x0 = xhit[zero]
            x1 = xhit[one]
            x2 = xhit[two]
            y0 = yhit[zero]
            y1 = yhit[one]
            y2 = yhit[two]
            a = (x1 - pp * x0) / (1 - pp)
            b = (y1 - pp * y0) / (1 - pp)
            c = (x2 - qq * x0) / (1 - qq)
            d = (y2 - qq * y0) / (1 - qq)
            rsquare = pp * ((x1 - x0) ** 2 + (y1 - y0) ** 2) / ((1 - pp) ** 2)
            ssquare = qq * ((x2 - x0) ** 2 + (y2 - y0) ** 2) / ((1 - qq) ** 2)
            e = c - a
            f = d - b
            if d == b:
                f = 0.000000001
            g = sqrt(e * e + f * f)
            k = 0.5 * (g * g + rsquare - ssquare) / g
            linelist0.append(-e / f)                            # coefficient of radical axis
            linelist1.append((a * e + b * f + g * k) / f)       # coefficient of radical axis

        linx, liny = CenterMassAlgorithm.reconstruct_common(p, x, y, z,
                                                            initial)  # in case the radical axis method does not deliver possible core positions
        subsets = itertools.combinations(statindex, 2)          # select pairs of radical axes

        xpointlist = []
        ypointlist = []
        for zero, one in subsets:
            a = linelist0[zero]
            b = linelist1[zero]
            c = linelist0[one]
            d = linelist1[one]
            aminc = a - c
            if a == c:
                aminc = 0.000000001             # to avoid singularity
            xint = (d - b) / aminc              # x coordinate of intersection of pair of radical axes
            yint = (a * d - b * c) / aminc      # y coordinate of intersection of pair of radical axes
            if abs(xint) < 600. and abs(yint)<600.:     # accept intersection point if not to far away
                xpointlist.append(xint)
                ypointlist.append(yint)

        if len(xpointlist) > 0:
            linx = mean(xpointlist)             # determine the average of the set of intersections of radical axes
            liny = mean(ypointlist)

        return linx, liny



class EllipsLdfAlgorithm(object):

    """ Special core and size estimator

    Estimates the core and shower size (electrons + muons) by a brute force
    inspection in a limited region, in combination with regression, with
    elliptic lateral densities in the neighborhood of cores found with both
    the CenterMassAlgorithm and AverageIntersectionAlgorithm.

    """

    @classmethod
    def reconstruct_common(cls, p, x, y, z=None, initial={}):
        """Reconstruct core position and shower size

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param z: height of detectors is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions: zenith and azimuth.

        """
        theta = initial.get('theta', 0.)
        phi = initial.get('phi', 0.)
        return cls.reconstruct(p, x, y, theta, phi)

    @classmethod
    def reconstruct(cls, p, x, y, theta, phi):
        """Reconstruct core position that performs best.
        Reconstruct shower size (electrons + muons

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param theta,phi: zenith and azimuth angle in rad.

        """

        xcmass, ycmass = CenterMassAlgorithm.reconstruct_common(p, x, y)
        chi2best = 10 ** 99
        factorbest = 1.
        gridsize = 20.          # determine within a grid around the barycenter estimation the position with minimum chi square value

        xbest, ybest, chi2best, factorbest = cls.selectbest(
            p, x, y, xcmass, ycmass, factorbest, chi2best, gridsize, theta, phi)

        gridsize = 5.            # determine within a smaller grid around the previous estimation the position with minimum chi square
        xbest1, ybest1, chi2best1, factorbest1 = cls.selectbest(
            p, x, y, xbest, ybest, factorbest, chi2best, gridsize, theta, phi)

        xlines, ylines = AverageIntersectionAlgorithm.reconstruct_common(p, x, y)
        chi2best = 10 ** 99
        factorbest = 1.
        gridsize = 50.          # determine within a grid around the radical axes estimation the position with minimum chi square

        xbest, ybest, chi2best, factorbest = cls.selectbest(
            p, x, y, xlines, ylines, factorbest, chi2best, gridsize, theta, phi)

        gridsize = 10.          # determine within a smaller grid around the previous estimation the position with minimum chi square
        xbest2, ybest2, chi2best2, factorbest2 = cls.selectbest(
            p, x, y, xbest, ybest, factorbest, chi2best, gridsize, theta, phi)


        if chi2best1 < chi2best2:       # select best estimation out of the two refined grid
            xbest, ybest, chi2best, factorbest = xbest1, ybest1, chi2best1, factorbest1
        else:
            xbest, ybest, chi2best, factorbest = xbest2, ybest2, chi2best2, factorbest2

        gridsize = 4.            # determine within a grid around the best estimation the position with minimum chi square
        core_x, core_y, chi2best, factorbest = cls.selectbest(
            p, x, y, xbest, ybest, factorbest, chi2best, gridsize, theta, phi)

        size = factorbest * ldf.EllipsLdf._Ne       # estimated shower size : number of electrons and muons !!!
        coefa = 0.524 * cos(theta) + 0.681
        coefb = 7.844 + 5.300 * cos(theta)
        enerpow = (log10(size) +coefb) / coefa      # estimated energy on the basis of relation between shower size and energy
        energy = 10 ** enerpow

        return core_x, core_y, chi2best, size, energy

    @staticmethod
    def selectbest(p, x, y, xstart, ystart, factorbest, chi2best, gridsize,
                   theta, phi):
        """selects the best core position in grid around (xstart, ystart) by
            optimizing the shower size (electr + muon) by means of regression.

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param xcmass,ycmass: start position of core in m.

        """
        xbest = xstart
        ybest = ystart

        a = ldf.EllipsLdf(zenith=theta, azimuth=phi)
        gridparam = 4
        for i in range(2 * gridparam + 1):
            xtry = xstart + (i - gridparam) * gridsize
            for j in range(2 * gridparam + 1):
                ytry = ystart + (j - gridparam) * gridsize
                xstations = array(x)
                ystations = array(y)
                r, angle = a.calculate_core_distance_and_angle(
                    xstations, ystations, xtry, ytry)
                rho = a.calculate_ldf_value(r, angle)

                mmdivl = 0.
                m = 0.
                l = 0.

                for ki, kj in zip(p, rho):
                    mmdivl += 1. * ki * ki / kj
                    m += ki
                    l += kj

                sizefactor = sqrt(mmdivl / l)
                chi2 = 2. * (sizefactor * l - m)

                if chi2 < chi2best:
                    factorbest = sizefactor
                    xbest = xtry
                    ybest = ytry
                    chi2best = chi2

        return xbest, ybest, chi2best, factorbest


class BruteForceAlgorithm(object):

    """ Brute force core and shower size (electrons + muons) estimator

    Estimates the core and shower size (electrons + muons) by a brute force
    inspection, in combination with regression, with elliptic lateral
    densities in the neighborhood of cores found with both the
    CenterMassAlgorithm and AverageIntersectionAlgorithm.

    NOT RECOMMENDED TO USE since it is extremely slow

    """

    @classmethod
    def reconstruct_common(cls, p, x, y, z=None, initial={}):
        """Reconstruct core position and shower size

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param z: height of detectors is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions: zenith and azimuth.

        """
        theta = initial.get('theta', 0.)
        phi = initial.get('phi', 0.)
        return cls.reconstruct(p, x, y, theta, phi)

    @classmethod
    def reconstruct(cls, p, x, y, theta, phi):
        """Reconstruct core position that performs best.
        Reconstruct showers size (electrons + muons) that fits best.

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param theta, phi: zenith and azimuth angle in rad.

        """
        chi2best = 10 ** 99
        factorbest = 1.
        xbest = 0.
        ybest = 0.
        gridsize = 10.            ################

        core_x, core_y, chi2best, factorbest = cls.selectbest(
            p, x, y, xbest, ybest, factorbest, chi2best, gridsize, theta, phi)

        size = factorbest * ldf.EllipsLdf._Ne
        coefa = 0.519 * cos(theta) + 0.684
        coefb = 7.84 +5.30 * cos(theta)
        enerpow = (log10(size) +coefb) / coefa
        energy = 10 ** enerpow

        return core_x, core_y, chi2best, size, energy

    @staticmethod
    def selectbest(p, x, y, xstart, ystart, factorbest, chi2best, gridsize,
                   theta, phi):
        """selects the best core position in grid around (xstart, ystart).

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param xcmass,ycmass: start position of core in m.

        """
        xbest = xstart
        ybest = ystart

        a = ldf.EllipsLdf(zenith=theta, azimuth=phi)
        gridparam = 75
        for i in range(2 * gridparam + 1):
            xtry = xstart + (i - gridparam) * gridsize
            for j in range(2 * gridparam + 1):
                ytry = ystart + (j - gridparam) * gridsize
                xstations = array(x)
                ystations = array(y)
                r, angle = a.calculate_core_distance_and_angle(
                    xstations, ystations, xtry, ytry)
                rho = a.calculate_ldf_value(r, angle)

                mmdivl = 0.
                m = 0.
                l = 0.

                for ki, kj in zip(p, rho):
                    mmdivl += 1. * ki * ki / kj
                    m += ki
                    l += kj

                sizefactor = sqrt(mmdivl / l)
                chi2 = 2. * (sizefactor * l - m)

                if chi2 < chi2best:
                    factorbest = sizefactor
                    xbest = xtry
                    ybest = ytry
                    chi2best = chi2

        return xbest, ybest, chi2best, factorbest

