""" PyTables table descriptions for data storage

    This module contains the table descriptions used by the detector
    simulation to store intermediate and final data in a HDF5 file.

"""
import tables


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


class EventObservables(tables.IsDescription):

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
    timestamp = tables.Time32Col()
    nanoseconds = tables.UInt32Col()
    ext_timestamp = tables.UInt64Col()

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


class Coincidence(tables.IsDescription):

    """Store information about a coincidence of stations within a cluster.

    An extensive air shower can trigger multiple stations, resulting in a set
    of events which are from the same shower.  This is called a coincidence.

    This table assigns an :attr:`id` to a coincidence and provides some
    additional information.  The events making up the coincidence can be looked
    up using the :class:`c_index` table.  Let ``coincidence`` be a row from
    this table, then you can do::

        >>> coincidence_id = coincidence['id']
        >>> event_ids = c_index[coincidence_id]
        >>> coincidence_event_list = [events[u] for u in event_ids]

    Note that all events included in the coincidence are not required to
    actually have measured particles.  For example, simulations include all
    events from the same shower in the coincidence, regardless of observed
    particles.  On the other hand, experimental datasets only include stations
    which have triggered, but may include events which have not actually
    measured the same shower, but simply measured other particles at the same
    time, by chance.

    Simulations may set the :attr:`r`, :attr:`phi`, :attr:`x`, :attr:`y` and
    :attr:`alpha` attributes to simulation parameters, like core position and
    cluster rotation.

    .. attribute:: id

        a unique identifier for the coincidence (only unique in this table)

    .. attribute:: N

        the number of triggered stations

    .. attribute:: r, phi, x, y

        The coordinates of the shower core in a simulation.

    .. attribute:: shower_theta, shower_phi

        The direction of the (simulated) shower.

    """
    id = tables.UInt32Col(pos=0)
    timestamp = tables.Time32Col(pos=1)
    nanoseconds = tables.UInt32Col(pos=2)
    ext_timestamp = tables.UInt64Col(pos=3)

    N = tables.UInt8Col(pos=4)
    x = tables.Float32Col(pos=5)
    y = tables.Float32Col(pos=6)
    zenith = tables.Float32Col(pos=7)
    azimuth = tables.Float32Col(pos=8)
    size = tables.Float32Col(pos=9)
    energy = tables.Float32Col(pos=10)


class ReconstructedEvent(tables.IsDescription):

    """Store information about reconstructed events"""

    # r, phi is core position

    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    alpha = tables.Float32Col()
    t1 = tables.Float32Col()
    t2 = tables.Float32Col()
    t3 = tables.Float32Col()
    t4 = tables.Float32Col()
    n1 = tables.Float32Col()
    n2 = tables.Float32Col()
    n3 = tables.Float32Col()
    n4 = tables.Float32Col()
    reference_theta = tables.Float32Col()
    reference_phi = tables.Float32Col()
    reconstructed_theta = tables.Float32Col()
    reconstructed_phi = tables.Float32Col()
    reference_core_pos = tables.Float32Col(shape=2)
    reconstructed_core_pos = tables.Float32Col(shape=2)
    reference_shower_size = tables.Float32Col()
    reconstructed_shower_size = tables.Float32Col()
    min_n134 = tables.Float32Col()


class KascadeEvent(tables.IsDescription):

    """Store events from KASCADE"""

    run_id = tables.IntCol()
    event_id = tables.Int64Col()
    timestamp = tables.Time32Col()
    nanoseconds = tables.UInt32Col()
    ext_timestamp = tables.UInt64Col()
    energy = tables.FloatCol()
    core_pos = tables.FloatCol(shape=2)
    zenith = tables.FloatCol()
    azimuth = tables.FloatCol()
    Num_e = tables.FloatCol()
    Num_mu = tables.FloatCol()
    dens_e = tables.FloatCol(shape=4)
    dens_mu = tables.FloatCol(shape=4)
    P200 = tables.FloatCol()
    T200 = tables.FloatCol()


class ReconstructedKascadeEvent(tables.IsDescription):

    """Store information about reconstructed events"""

    # r, phi is core position

    id = tables.UInt32Col()
    station_id = tables.UInt8Col()
    r = tables.Float32Col()
    phi = tables.Float32Col()
    alpha = tables.Float32Col()
    t1 = tables.Float32Col()
    t2 = tables.Float32Col()
    t3 = tables.Float32Col()
    t4 = tables.Float32Col()
    n1 = tables.Float32Col()
    n2 = tables.Float32Col()
    n3 = tables.Float32Col()
    n4 = tables.Float32Col()
    reference_theta = tables.Float32Col()
    reference_phi = tables.Float32Col()
    reconstructed_theta = tables.Float32Col()
    reconstructed_phi = tables.Float32Col()
    min_n134 = tables.Float32Col()

    k_energy = tables.FloatCol()
    k_core_pos = tables.FloatCol(shape=2)
    k_Num_e = tables.FloatCol()
    k_Num_mu = tables.FloatCol()
    k_dens_e = tables.FloatCol(shape=4)
    k_dens_mu = tables.FloatCol(shape=4)
    k_P200 = tables.FloatCol()
    k_T200 = tables.FloatCol()
