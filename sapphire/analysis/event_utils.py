""" Get data from HiSPARC events

This module contains functions to derive data from HiSPARC events.
Common tasks for data reconstruction are getting the particle density
and shower arrival time in detectors or a station. These functions are
aware of processed events (i.e. reconstructed number of MIPs, arrival
times and trigger time) and stations.

"""
from ..utils import ERR

from numpy import nan, nanmin, nanmean


NO_OFFSET = [0., 0., 0., 0.]


def station_density(event, detector_ids=None, station=None):
    """Get particle density in a station

    Detectors with error values will be ignored.

    :param event: Processed event row.
    :param detector_ids: list of detectors ids to consider. If None, the
        detectors in the station object will be used.
    :param station: :class:`~sapphire.clusters.Station` object.
    :return: average density over the chosen detectors.

    """
    if detector_ids is None:
        detector_ids = get_detector_ids(station, event)
    p = nanmean(detector_densities(event, detector_ids=detector_ids,
                                   station=station))
    return p


def detector_densities(event, detector_ids=None, station=None):
    """Get particle density in station detectors

    :param event: Processed event row.
    :param detector_ids: list of detectors ids for which to get particle
        densities.
    :param station: :class:`~sapphire.clusters.Station` object.
    :return: density in each chosen detector.

    """
    if detector_ids is None:
        detector_ids = get_detector_ids(station, event)
    p = [detector_density(event, id, station) for id in detector_ids]
    return p


def detector_density(event, detector_id, station=None):
    """Get particle density in station detector

    :param event: Processed event row.
    :param detector_id: detector id for which to get particle density.
    :param station: :class:`~sapphire.clusters.Station` object, used to
        determine the detector size.
    :return: density in the chosen detector.

    """
    number_of_particles = event['n%d' % (detector_id + 1)]
    try:
        area = station.detectors[detector_id].get_area()
    except (AttributeError, IndexError):
        area = .5
    if number_of_particles not in ERR:
        p = number_of_particles / area
    else:
        p = nan
    return p


def station_arrival_time(event, reference_ext_timestamp,
                         detector_ids=None, offsets=NO_OFFSET, station=None):
    """Get station arrival time, i.e. first detector hit

    Arrival time of first detector hit in the station. The returned time
    is relative to reference_ext_timestamp, because floats do not have
    enough precision to represent large timestamps in nanoseconds.

    :param event: Processed event row.
    :param reference_ext_timestamp: reference extended timestamp (in ns).
        The returned station arrival time will be relative to this timestamp.
        Often best to use the timestamp of the first event in a coincidence.
    :param detector_ids: list of detectors ids for which to consider.
    :param offsets: list of detector time offsets.
    :param station: :class:`~sapphire.clusters.Station` object, used to
        determine the number of detectors.
    :return: shower arrival time of the station relative to the
             reference timestamp.

    """
    if detector_ids is None:
        detector_ids = get_detector_ids(station, event)
    if event['t_trigger'] in ERR:
        t = nan
    else:
        t_first = nanmin(detector_arrival_times(event, detector_ids, offsets,
                                                station))
        t = ((int(event['ext_timestamp']) - int(reference_ext_timestamp)) -
             event['t_trigger'] + t_first)
    return t


def relative_detector_arrival_times(event, reference_ext_timestamp,
                                    detector_ids=None, offsets=NO_OFFSET,
                                    station=None):
    """Get relative arrival times for all detectors

    :param event: Processed event row.
    :param reference_ext_timestamp: reference extended timestamp (in ns).
        The returned station arrival time will be relative to this timestamp.
        Often best to use the timestamp of the first event in a coincidence.
    :param detector_ids: list of detectors ids for which to get arrival times.
    :param offsets: list of detector time offsets.
    :param station: :class:`~sapphire.clusters.Station` object, used to
        determine the number of detectors.
    :return: list of shower arrival times relative to the given reference.

    """
    if detector_ids is None:
        detector_ids = get_detector_ids(station, event)
    if event['t_trigger'] in ERR:
        t = [nan] * len(detector_ids)
    else:
        arrival_times = detector_arrival_times(event, detector_ids,
                                               offsets, station)
        t = [(int(event['ext_timestamp']) - int(reference_ext_timestamp)) -
             event['t_trigger'] + arrival_time
             for arrival_time in arrival_times]
    return t


def detector_arrival_times(event, detector_ids=None, offsets=NO_OFFSET,
                           station=None):
    """Get corrected arrival times for all detectors

    :param event: Processed event row.
    :param detector_ids: list of detectors ids for which to get arrival times.
    :param offsets: list of detector time offsets.
    :param station: :class:`~sapphire.clusters.Station` object, used to
        determine the number of detectors.
    :return: list of shower arrival times relative to the start of the trace.

    """
    if detector_ids is None:
        detector_ids = get_detector_ids(station, event)
    t = [detector_arrival_time(event, id, offsets) for id in detector_ids]
    return t


def detector_arrival_time(event, detector_id, offsets=NO_OFFSET):
    """Get corrected arrival time for a detector

    :param event: Event row.
    :param detector_id: detector id for which to get arrival times.
    :param offsets: list of detector time offsets.
    :return: arrival time corrected by the offset.

    """
    arrival_time = event['t%d' % (detector_id + 1)]
    if arrival_time not in ERR:
        t = arrival_time - offsets[detector_id]
    else:
        t = nan
    return t


def get_detector_ids(station=None, event=None):
    """Determine the detector ids based on the station object or event data

    Returns a list of detectors that should be present.

    Note: Event based determination might not work for simulated data,
    since it currently does not simulate pulseheights.

    :param event: Event row.
    :param station: :class:`~sapphire.clusters.Station` object.
    :return: list of detector_ids.

    """
    if station is not None:
        detector_ids = range(len(station.detectors))
    elif event is not None:
        detector_ids = [i for i, ph in enumerate(event['pulseheights'])
                        if ph != -1]
    else:
        detector_ids = range(4)
    return detector_ids
