""" PyTables table descriptions for data storage

    This module contains the table descriptions used by the detector
    simulation to store intermediate and final data in a HDF5 file.

"""
import tables


class Particle(tables.IsDescription):
    """Store information about shower particles reaching round level"""
    id = tables.UInt32Col()
    pid = tables.Int8Col()
    core_distance = tables.Float32Col()
    polar_angle = tables.Float32Col()
    x = tables.Float32Col()
    y = tables.Float32Col()
    arrival_time = tables.Float32Col()
    energy = tables.Float32Col()


class SimulationHeader(tables.IsDescription):
    """Header storing single station information of a simulated event"""
    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    alpha = tables.Float32Col()


class ParticleEvent(tables.IsDescription):
    """Store information about the particles hitting a detector"""
    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    detector_id = tables.UInt8Col()
    pid = tables.Int8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    time = tables.Float32Col()
    energy = tables.Float32Col()


class ObservableEvent(tables.IsDescription):
    """Store information about the observables of an event"""
    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    x = tables.Float32Col()
    y = tables.Float32Col()
    alpha = tables.Float32Col()
    N = tables.UInt8Col()
    t1 = tables.Float32Col()
    t2 = tables.Float32Col()
    t3 = tables.Float32Col()
    t4 = tables.Float32Col()
    n1 = tables.UInt16Col()
    n2 = tables.UInt16Col()
    n3 = tables.UInt16Col()
    n4 = tables.UInt16Col()


class CoincidenceEvent(tables.IsDescription):
    """Store information about a coincidence"""
    id = tables.UInt32Col()
    N = tables.UInt8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    x = tables.Float32Col()
    y = tables.Float32Col()
    alpha = tables.Float32Col()
