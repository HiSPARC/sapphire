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
from itertools import izip, izip_longest, combinations
import warnings

from numpy import isnan, nan, cos, sqrt, log10, mean, array, linspace

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

    def reconstruct_event(self, event, detector_ids=None, initial={}):
        """Reconstruct a single event

        :param event: an event (e.g. from an events table), or any
            dictionary-like object containing the keys necessary for
            reconstructing the direction of a shower (e.g. number of mips).
        :param detector_ids: list of the detectors to use for
            reconstruction. The detector ids are 0-based, unlike the
            column names in the esd data.
        :param initial: dictionary with already reconstructed shower
                        parameters.
        :return: (x, y) core position in m, shower size, and energy in eV.

        """
        p, x, y, z = ([], [], [], [])
        if detector_ids is None:
            detector_ids = range(4)
        self.station.cluster.set_timestamp(event['timestamp'])
        for id in detector_ids:
            p_detector = detector_density(event, id, self.station)
            if not isnan(p_detector):
                dx, dy, dz = self.station.detectors[id].get_coordinates()
                p.append(p_detector)
                x.append(dx)
                y.append(dy)
                z.append(dz)
        if len(p) >= 3:
            core_x, core_y, size, energy = \
                self.estimator.reconstruct_common(p, x, y, z, initial)
        else:
            core_x, core_y, size, energy = (nan, nan, nan, nan)
        return core_x, core_y, size, energy

    def reconstruct_events(self, events, detector_ids=None, progress=True,
                           initials=[]):
        """Reconstruct events

        :param events: the events table for the station from an ESD data
                       file.
        :param detector_ids: detectors which use for the reconstructions.
        :param progress: if True shows a progress bar.
        :param initials: list of dictionaries with already reconstructed shower
                         parameters.
        :return: (x, y) core positions in m, shower sizes, and energies in eV.

        """
        events = pbar(events, show=progress)
        events_init = izip_longest(events, initials)
        cores = [self.reconstruct_event(event, detector_ids, initial)
                 for event, initial in events_init]
        if len(cores):
            core_x, core_y, size, energy = zip(*cores)
        else:
            core_x, core_y, size, energy = ((), (), (), ())
        return core_x, core_y, size, energy


class CoincidenceCoreReconstruction(object):

    """Reconstruct core for coincidences

    This class is aware of 'coincidences' and 'clusters'.  Initialize
    this class with a 'cluster' and you can reconstruct a coincidence
    using :meth:`reconstruct_coincidence`.

    :param cluster: :class:`sapphire.clusters.BaseCluster` object.

    """

    def __init__(self, cluster):
        self.estimator = EllipseLdfAlgorithm
        self.cluster = cluster

    def reconstruct_coincidence(self, coincidence, station_numbers=None,
                                initial={}):
        """Reconstruct a single coincidence

        :param coincidence: a coincidence list consisting of
                            multiple (station_number, event) tuples
        :param station_numbers: list of station numbers, to only use
                                events from those stations.
        :param initial: dictionary with already reconstructed shower
                        parameters.
        :param initial: dictionary with initial data.
        :return: (x, y) core position in m, shower size, and energy in eV.

        """
        p, x, y, z = ([], [], [], [])

        try:
            self.cluster.set_timestamp(coincidence[0][1]['timestamp'])
        except IndexError:
            return (nan, nan)

        for station_number, event in coincidence:
            if station_numbers is not None:
                if station_number not in station_numbers:
                    continue
            station = self.cluster.get_station(station_number)
            p_station = station_density(event, range(4), station)
            if not isnan(p_station):
                sx, sy, sz = station.calc_center_of_mass_coordinates()
                p.append(p_station)
                x.append(sx)
                y.append(sy)
                z.append(sz)

        if len(p) >= 3:
            core_x, core_y, size, energy = \
                self.estimator.reconstruct_common(p, x, y, z, initial)
        else:
            core_x, core_y, size, energy = (nan, nan, nan, nan)
        return core_x, core_y, size, energy

    def reconstruct_coincidences(self, coincidences, station_numbers=None,
                                 progress=True, initials=[]):
        """Reconstruct all coincidences

        :param coincidences: a list of coincidences, each consisting of
                             multiple (station_number, event) tuples.
        :param station_numbers: list of station numbers, to only use
                                events from those stations.
        :param progress: if True shows a progress bar.
        :param initials: list of dictionaries with already reconstructed shower
                         parameters.
        :return: (x, y) core positions in m, shower sizes, and energies in eV.

        """
        coincidences = pbar(coincidences, show=progress)
        coin_init = izip_longest(coincidences, initials)
        cores = [self.reconstruct_coincidence(coincidence, station_numbers,
                                              initial)
                 for coincidence, initial in coin_init]
        if len(cores):
            core_x, core_y, size, energy = zip(*cores)
        else:
            core_x, core_y, size, energy = ((), (), (), ())
        return core_x, core_y, size, energy


class CoincidenceCoreReconstructionDetectors(
        CoincidenceCoreReconstruction):

    """Reconstruct core for coincidences using each detector

    Instead of using the average station particle density this class
    uses the particle density in each detector for the reconstruction.

    """

    def reconstruct_coincidence(self, coincidence, station_numbers=None,
                                initial={}):
        """Reconstruct a single coincidence

        :param coincidence: a coincidence list consisting of
                            multiple (station_number, event) tuples
        :param station_numbers: list of station numbers, to only use
                                events from those stations.
        :param initial: dictionary with already reconstructed shower
                        parameters.
        :return: (x, y) core position in m.

        """
        p, x, y, z = ([], [], [], [])

        try:
            self.cluster.set_timestamp(coincidence[0][1]['timestamp'])
        except IndexError:
            return (nan, nan)

        for station_number, event in coincidence:
            if station_numbers is not None:
                if station_number not in station_numbers:
                    continue
            station = self.cluster.get_station(station_number)
            for id in range(4):
                p_detector = detector_density(event, id, station)
                if not isnan(p_detector):
                    dx, dy, dz = station.detectors[id].get_coordinates()
                    p.append(p_detector)
                    x.append(dx)
                    y.append(dy)
                    z.append(dz)

        if len(p) >= 3:
            core_x, core_y = self.estimator.reconstruct_common(p, x, y, z,
                                                               initial)
        else:
            core_x, core_y = (nan, nan)
        return core_x, core_y


class CenterMassAlgorithm(object):

    """Simple core estimator

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
        :return: (x, y) core position in m. Shower size and energy are nan.

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
        :return: (x, y) core position in m. Shower size and energy are nan.

        """
        core_x = sum(density * xi for density, xi in zip(p, x)) / sum(p)
        core_y = sum(density * yi for density, yi in zip(p, y)) / sum(p)
        return core_x, core_y, nan, nan


class AverageIntersectionAlgorithm(object):

    """Core estimator

    To the densities in 3 stations correspond 2 possible cores. The line
    through these points is quite stable for the lateral distribution function.
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
        :return: (x, y) core position in m. Shower size and energy are nan.

        """
        if len(p) < 4 or len(x) < 4 or len(y) < 4:
            warnings.warn('This algorithm requires at least 4 detections.',
                          UserWarning)
            return nan, nan, nan, nan

        theta = initial.get('theta', nan)
        if not isnan(theta):
            p = [density * cos(theta) for density in p]

        return cls.reconstruct(p, x, y)

    @classmethod
    def reconstruct(cls, p, x, y):
        """Calculate center of mass

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :return: (x, y) core position in m. Shower size and energy are nan.

        """
        phit = []
        xhit = []
        yhit = []
        for i in range(len(p)):
            if p[i] > 0.001:
                phit.append(p[i])
                xhit.append(x[i])
                yhit.append(y[i])

        statindex = range(len(phit))

        linelist0 = []
        linelist1 = []
        # select triples of stations
        for zero, one, two in combinations(statindex, 3):
            a, b, rsquare = cls.calculate(phit, xhit, yhit, zero, one)
            c, d, ssquare = cls.calculate(phit, xhit, yhit, zero, two)
            e = c - a
            if d == b:
                f = 0.000000001
            else:
                f = d - b
            g = sqrt(e ** 2 + f ** 2)
            k = 0.5 * (g * g + rsquare - ssquare) / g
            # coefficients of radical axis
            linelist0.append(-e / f)
            linelist1.append((a * e + b * f + g * k) / f)

        xpointlist = []
        ypointlist = []
        # select pairs of radical axes
        for zero, one in combinations(statindex, 2):
            a = linelist0[zero]
            b = linelist1[zero]
            c = linelist0[one]
            d = linelist1[one]
            if a == c:
                aminc = 0.000000001
            else:
                aminc = a - c
            # x and y coordinates of intersection of pair of radical axes
            xint = (d - b) / aminc
            yint = (a * d - b * c) / aminc
            # accept intersection point if not to far away
            if abs(xint) < 600. and abs(yint) < 600.:
                xpointlist.append(xint)
                ypointlist.append(yint)

        if len(xpointlist):
            # determine the average of the set of intersections of radical axes
            core_x = mean(xpointlist)
            core_y = mean(ypointlist)
        else:
            # Fallback to CenterMass
            return CenterMassAlgorithm.reconstruct_common(p, x, y)

        return core_x, core_y, nan, nan

    @staticmethod
    def calculate(p, x, y, i, j):
        """Perform a calculation that is used multiple times"""

        m = 2.3  # optimized value in powerlaw  r ^(-m)  for density

        pp = (p[i] / p[j]) ** (2. / m)
        if pp == 1:
            pp = 1.000001
        a = (x[j] - pp * x[i]) / (1 - pp)
        b = (y[j] - pp * y[i]) / (1 - pp)
        square = (pp * ((x[j] - x[i]) ** 2 + (y[j] - y[i]) ** 2) /
                  ((1 - pp) ** 2))
        return a, b, square


class EllipseLdfAlgorithm(object):

    """Core and size estimator using an LDF

    Estimates the core and shower size (electrons + muons) by a brute force
    inspection in a limited region, in combination with regression, with
    elliptic lateral densities in the neighbourhood of cores found with both
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
        :return: (x, y) core position in m, shower size, and energy.

        """
        theta = initial.get('theta', 0.)
        phi = initial.get('phi', 0.)
        return cls.reconstruct(p, x, y, theta, phi)

    @classmethod
    def reconstruct(cls, p, x, y, theta, phi):
        """Reconstruct best fitting shower core position, size, and energy

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param theta,phi: zenith and azimuth angle in rad.
        :return: (x, y) core position in m, shower size, and energy.

        """
        xcm, ycm, _, _ = CenterMassAlgorithm.reconstruct_common(p, x, y)
        chi2cm = 10 ** 99
        factorcm = 1.

        for gridsize in [20., 5.]:
            xcm, ycm, chi2cm, factorcm = cls.selectbest(
                p, x, y, xcm, ycm, factorcm, chi2cm, gridsize, theta, phi)

        xai, yai, _, _ = \
            AverageIntersectionAlgorithm.reconstruct_common(p, x, y)
        chi2ai = 10 ** 99
        factorai = 1.

        for gridsize in [50., 10.]:
            xai, yai, chi2ai, factorai = cls.selectbest(
                p, x, y, xai, yai, factorai, chi2ai, gridsize, theta, phi)

        if chi2cm < chi2ai:
            xbest, ybest, chi2best, factorbest = xcm, ycm, chi2cm, factorcm
        else:
            xbest, ybest, chi2best, factorbest = xai, yai, chi2ai, factorai

        # determine within a grid around the best estimation the position
        # with minimum chi square
        gridsize = 4.
        core_x, core_y, chi2best, factorbest = cls.selectbest(
            p, x, y, xbest, ybest, factorbest, chi2best, gridsize, theta, phi)
        # estimated shower size, number of electrons and muons
        size = factorbest * ldf.EllipseLdf._Ne
        coefa = 0.524 * cos(theta) + 0.681
        coefb = 7.844 + 5.30 * cos(theta)
        # estimated energy based on relation between shower size and energy
        enerpow = (log10(size) + coefb) / coefa
        energy = 10 ** enerpow

        return core_x, core_y, size, energy

    @staticmethod
    def selectbest(p, x, y, xstart, ystart, factorbest, chi2best, gridsize,
                   theta, phi):
        """Select best core position in grid around (xstart, ystart)

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param xstart,ystart: start position of core in m.
        :param factorbest: current best estimate for shower scale (x ldf.Ne).
        :param chi2best: chi2 for current core and size.
        :param gridsize: size of grid around current core.
        :param theta,phi: shower direction.
        :return: core position, chi2 and shower scale.

        """
        xstations = array(x)
        ystations = array(y)

        xbest = xstart
        ybest = ystart

        elldf = ldf.EllipseLdf(zenith=theta, azimuth=phi)
        gridparam = 4
        gridedge = gridparam * gridsize
        gridpoints = linspace(-gridedge, gridedge, gridparam * 2 + 1)
        for xtry in gridpoints + xstart:
            for ytry in gridpoints + ystart:
                r, angle = elldf.calculate_core_distance_and_angle(
                    xstations, ystations, xtry, ytry)
                rho = elldf.calculate_ldf_value(r, angle)

                mmdivl = 0.
                m = 0.
                l = 0.

                for ki, kj in zip(p, rho):
                    mmdivl += 1. * ki * ki / kj
                    m += ki
                    l += kj

                sizefactor = sqrt(mmdivl / l)
                with warnings.catch_warnings(record=True):
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

    .. warning::
        NOT RECOMMENDED TO USE since it is extremely slow

    """

    @classmethod
    def reconstruct_common(cls, p, x, y, z=None, initial={}):
        """Reconstruct core position and shower size

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param z: height of detectors is ignored.
        :param initial: dictionary containing values from previous
                        reconstructions, i.e. zenith (theta) and azimuth (phi).

        """
        theta = initial.get('theta', 0.)
        phi = initial.get('phi', 0.)
        return cls.reconstruct(p, x, y, theta, phi)

    @classmethod
    def reconstruct(cls, p, x, y, theta, phi):
        """Reconstruct best fitting shower core position, size, and energy

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param theta,phi: zenith and azimuth angle in rad.
        :return: (x, y) core position in m, shower size, and energy.

        """
        chi2best = 10 ** 99
        factorbest = 1.
        xbest = 0.
        ybest = 0.
        gridsize = 10.

        core_x, core_y, chi2best, factorbest = cls.selectbest(
            p, x, y, xbest, ybest, factorbest, chi2best, gridsize, theta, phi)

        size = factorbest * ldf.EllipseLdf._Ne
        coefa = 0.519 * cos(theta) + 0.684
        coefb = 7.84 + 5.30 * cos(theta)
        enerpow = (log10(size) + coefb) / coefa
        energy = 10 ** enerpow

        return core_x, core_y, size, energy

    @staticmethod
    def selectbest(p, x, y, xstart, ystart, factorbest, chi2best, gridsize,
                   theta, phi):
        """Select the best core position in grid around (xstart, ystart).

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param xstart,ystart: start position of core in m.
        :param factorbest: current best estimate for shower scale (x ldf.Ne).
        :param chi2best: chi2 for current core and size.
        :param gridsize: size of grid around current core.
        :param theta,phi: shower direction.
        :return: core position, chi2 and shower scale.

        """
        xstations = array(x)
        ystations = array(y)

        xbest = xstart
        ybest = ystart

        elldf = ldf.EllipseLdf(zenith=theta, azimuth=phi)
        gridparam = 75
        gridedge = gridparam * gridsize
        gridpoints = linspace(-gridedge, gridedge, gridparam * 2 + 1)
        for xtry in gridpoints + xstart:
            for ytry in gridpoints + ystart:
                r, angle = elldf.calculate_core_distance_and_angle(
                    xstations, ystations, xtry, ytry)
                rho = elldf.calculate_ldf_value(r, angle)

                mmdivl = 0.
                m = 0.
                l = 0.

                for ki, kj in zip(p, rho):
                    mmdivl += ki ** 2 / kj
                    m += ki
                    l += kj

                sizefactor = sqrt(mmdivl / l)
                with warnings.catch_warnings(record=True):
                    chi2 = 2. * (sizefactor * l - m)
                if chi2 < chi2best:
                    factorbest = sizefactor
                    xbest = xtry
                    ybest = ytry
                    chi2best = chi2

        return xbest, ybest, chi2best, factorbest
