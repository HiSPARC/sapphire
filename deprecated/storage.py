""" PyTables table descriptions for data storage

    This module contains the table descriptions used by the detector
    simulation to store intermediate and final data in a HDF5 file.

"""
import tables


class SimulationEventHeader(tables.IsDescription):
    """Header storing single station information of a simulated event

    Some simulations do not only write :class:`SimulationEventObservables`, but
    also raw simulation data.  This may be individual particles hitting
    detectors, for example.  To keep track of the simulation, headers may be
    written in this table.

    In the case of the groundparticles simulation, one header is written for
    each simulated event with an :attr:`station_id` set to 0, containing event
    information.  Additionally, for each station a header is written (with the
    :attr:`station_id` set to the station identifier) with station information.

    .. attribute:: id

        an event identifier (unique in this table)

    .. attribute:: station_id

        station identifier, such that you can do::

            >>> station = cluster.stations[station_id]

    .. attribute:: r, phi

        dataset-specific.  This might be the coordinates of the shower core in
        a simulation, or the location of the cluster center.  Consult the
        documentation provided by the simulation code.

    .. attribute:: alpha

        dataset-specific.  This might be the rotation of the cluster around
        its center.  Consult the documentation provided by the simulation code.

    """
    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    alpha = tables.Float32Col()


class SimulationParticle(tables.IsDescription):
    """Store information about the particles hitting a detector

    Simulations which track individual particles write particle information in
    this table.  Position, arrival time and energy, as well as the detector
    which detected this particle are stored.

    .. attribute:: id

        a unique identifier for the simulated event (only unique in this table)

    .. attribute:: station_id

        station identifier, such that you can do::

            >>> station = cluster.stations[station_id]

    .. attribute:: detector_id

        detector identifier, such that you can do::

            >>> station = cluster.stations[station_id]
            >>> detector = station.detectors[detector_id]

    .. attribute:: pid

        a particle identifier. Possible values are determined by the
        simulation package.

    .. attribute:: r, phi

        particle position in polar coordinates

    .. attribute:: time

        arrival time of the particle [ns]

    .. attribute:: energy

        particle energy [GeV]

    """
    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    detector_id = tables.UInt8Col()
    pid = tables.Int8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    time = tables.Float32Col()
    energy = tables.Float32Col()


class SimulationEventObservables(tables.IsDescription):

    """Store information about the observables of an event.

    The observables are described for each station independently.  So, for each
    event (with a unique :attr:`id`), there is a table row for each station
    (with a unique :attr:`station_id`), such that only the (id, station_id)
    combinations are unique in the table.

    .. attribute:: id

        a unique identifier for the simulated event (only unique in this table)

    .. attribute:: station_id

        station identifier, such that you can do::

            >>> station = cluster.stations[station_id]

    .. attribute:: r, phi, x, y

        coordinates of the station.  Depending on the simulation, this might be
        constant throughout the simulation, or it might change event by event.

    .. attribute:: alpha

        rotation of the station around its center

    .. attribute:: N

        number of detectors with at least one particle

    """
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
    n1 = tables.Float32Col()
    n2 = tables.Float32Col()
    n3 = tables.Float32Col()
    n4 = tables.Float32Col()


class ShowerParticle(tables.IsDescription):
    """Store information about shower particles reaching round level

    This table stores particles from shower simulations.  For example, AIRES
    simulations produce ``grdpcles`` files containing all particles which
    reached ground level.  These files can be read and their contents can be
    stored in this table.

    .. attribute:: id

        a unique identifier for the particle (unique in this table)

    .. attribute:: pid

        a particle identifier. Possible values are determined by the
        simulation package.

    .. attribute:: core_distance

        distance from the particle position to the shower core

    .. attribute:: polar_angle

        angle of the particle position vector to a reference line

    .. attribute:: x, y

        particle position

    .. attribute:: arrival_time

        arrival time of the particle [ns]

    .. attribute:: energy

        particle energy [GeV]

    """
    id = tables.UInt32Col()
    pid = tables.Int8Col()
    core_distance = tables.Float32Col()
    polar_angle = tables.Float32Col()
    x = tables.Float32Col()
    y = tables.Float32Col()
    arrival_time = tables.Float32Col()
    energy = tables.Float32Col()
