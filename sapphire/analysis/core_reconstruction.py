""" Core reconstruction

    This module contains two classes that can be used to reconstruct
    HiSPARC events and coincidences. The classes know how to extract the
    relevant information from the station and event or cluster and
    coincidence. Various algorithms which do the reconstruction are also
    defined here. The algorithms require positions and particle densties to
    do the reconstruction.

    Each algorithm has a :meth:`~CenterMassAlgorithm.reconstruct_common`
    method which always requires the same arguments: particle denisties,
    x, and y positions, optionally z positions can be given. The data is
    then prepared for the algorithm and passed to the `reconstruct`
    method which returns the reconstructed x and y coordinates.

"""
from __future__ import division

from numpy import nan

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
    def reconstruct_common(cls, p, x, y, z=None):
        """Reconstruct core position

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param z: height of detectors is ignored.

        """
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
