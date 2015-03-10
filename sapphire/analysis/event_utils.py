""" Get data from HiSPARC events

This module contains functions to derive data from HiSPARC events.

"""
from ..utils import ERR

from numpy import nan, nanmin


NO_OFFSET = [0., 0., 0., 0.]


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
    :param station: Station object, used to determine the numnber of detectors.
    :returns: list of shower arrival times relative to the given reference.

    """
    if detector_ids is None:
        detector_ids = range(len(station.detectors))
    if event['t_trigger'] in ERR:
        t = [nan] * len(detector_ids)
    else:
        arrival_times = detector_arrival_times(event, detector_ids,
                                               offsets, station)
        t = [(int(event['ext_timestamp']) - int(reference_ext_timestamp)) -
             event['t_trigger'] + arrival_time
             for arrival_time in arrival_times]
    return t


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
    :param station: Station object, used to determine the numnber of detectors.
    :returns: shower arrival time of the station relative to the
              reference timestamp.

    """
    if detector_ids is None:
        detector_ids = range(len(station.detectors))
    if event['t_trigger'] in ERR:
        t = nan
    else:
        t_first = nanmin(detector_arrival_times(event, detector_ids, offsets,
                                                station))
        t = ((int(event['ext_timestamp']) - int(reference_ext_timestamp)) -
             event['t_trigger'] + t_first)
    return t


def detector_arrival_times(event, detector_ids=None, offsets=NO_OFFSET,
                           station=None):
    """Get corrected arrival times for all detectors

    :param event: Processed event row.
    :param detector_ids: list of detectors ids for which to get arrival times.
    :param offsets: list of detector time offsets.
    :param station: Station object, used to determine the numnber of detectors.
    :returns: list of shower arrival times relative to the start of the trace.

    """
    if detector_ids is None:
        detector_ids = range(len(station.detectors))
    t = [detector_arrival_time(event, id, offsets) for id in detector_ids]
    return t


def detector_arrival_time(event, detector_id, offsets=NO_OFFSET):
    """Get corrected arrival time for a detector

    :param event: Event row.
    :param detector_id: detector id for which to get arrival times.
    :param offsets: list of detector time offsets.
    :returns: arrival time corrected by the offset.

    """
    arrival_time = event['t%d' % (detector_id + 1)]
    if arrival_time not in ERR:
        t = arrival_time - offsets[detector_id]
    else:
        t = nan
    return t
