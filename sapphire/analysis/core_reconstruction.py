from __future__ import division


class CenterMassAlgorithm(object):

    """ Simple core estimator

    Estimates the core by center of mass of the measurements.

    """

    @classmethod
    def reconstruct_common(cls, p, x, y, z=None):
        """Reconstruct core

        :param p: detector particle density in m^-2.
        :param x,y: positions of detectors in m.
        :param z: height of detectors is ignored.

        """
        core_x = sum(density * xi for density, xi in zip(p, x)) / sum(p)
        core_y = sum(density * yi for density, yi in zip(p, y)) / sum(p)
        return core_x, core_y


class EventReconstruction(CenterMassAlgorithm):

    """Reconstruct core for station events

    This class is aware of 'events' and 'stations'.  Initialize this class
    with a 'station' and you can reconstruct events using
    :meth:`reconstruct_event`.

    :param station: :class:`sapphire.clusters.Station` object.

    """

    def __init__(self, station):
        self.station = station
        detectors = [d.get_coordinates() for d in self.station.detectors]
        self.area = [d.get_area() for d in self.station.detectors]
        self.x, self.y, self.z = zip(*detectors)

    def reconstruct_event(self, event, detector_ids=[0, 1, 2, 3]):
        """Reconstruct a single event

        :param event: an event (e.g. from an events table), or any
            dictionary-like object containing the keys necessary for
            reconstructing the direction of a shower (e.g. number of mips).
        :param detector_ids: list of the detectors to use for
            reconstruction. The detector ids are 0-based, unlike the
            column names in the esd data.
        :returns: x, y core position in m.

        """
        p = [event['n%d' % (id + 1)] / self.area[id] for id in detector_ids]
        x = [self.x[id] for id in detector_ids]
        y = [self.y[id] for id in detector_ids]
        z = [self.z[id] for id in detector_ids]
        if len([i for i in n if i not in ERR]) >= 3:
            x, y = self.reconstruct_common(p, x, y, z)
        else:
            x = nan
            y = nan
        return x, y

    def reconstruct_events(self, events, detector_ids=[0, 1, 2, 3]):
        """Reconstruct events

        :param events: the events table for the station from an ESD data
                       file.
        :param detector_ids: detectors which use for the reconstructions.
        :returns: x, y core positions in m.

        """
        cores = [self.reconstruct_event(event, detector_ids)
                 for event in pbar(events)]
        x, y = zip(*cores)
        return x, y

