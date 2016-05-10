""" PyTables table descriptions for data storage

    This module contains the table descriptions used by the detector
    simulation to store intermediate and final data in a HDF5 file.

"""
import tables


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
    up using the ``c_index`` table.  Let ``coincidence`` be a row from
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

    Simulations may set the :attr:`x`, :attr:`y`, :attr:`zenith`,
    :attr:`azimuth`, :attr:`size`, and :attr:`energy` attributes to simulation
    parameters, like core position and shower parameters.

    .. attribute:: id

        a unique identifier for the coincidence (only unique in this table)

    .. attribute:: N

        the number of triggered stations

    .. attribute:: x

        The x coordinate of the shower core in a simulation.

    .. attribute:: y

        The y coordinate of the shower core in a simulation.

    .. attribute:: zenith

        The zenith direction of the (simulated) shower.

    .. attribute:: azimuth

        The azimuth direction of the (simulated) shower.

    .. attribute:: size

        The size (number of leptons) of the (simulated) shower.

    .. attribute:: energy

        The primary particle energy of the (simulated) shower.

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


class TimeDelta(tables.IsDescription):

    """Store time differences"""

    ext_timestamp = tables.UInt64Col(pos=0)
    timestamp = tables.UInt32Col(pos=1)
    nanoseconds = tables.UInt32Col(pos=2)
    delta = tables.FloatCol(pos=3)


class ReconstructedCoincidence(tables.IsDescription):

    """Store information about reconstructed coincidences"""

    id = tables.UInt32Col(pos=1)
    ext_timestamp = tables.UInt64Col(pos=2)
    min_n = tables.Float32Col(pos=3)

    x = tables.Float32Col(pos=4)
    y = tables.Float32Col(pos=5)
    zenith = tables.Float32Col(pos=6)
    azimuth = tables.Float32Col(pos=7)
    size = tables.Float32Col(pos=8)
    energy = tables.Float32Col(pos=9)
    error_x = tables.Float32Col(pos=10)
    error_y = tables.Float32Col(pos=11)
    error_zenith = tables.Float32Col(pos=12)
    error_azimuth = tables.Float32Col(pos=13)
    error_size = tables.Float32Col(pos=14)
    error_energy = tables.Float32Col(pos=15)

    reference_x = tables.Float32Col(pos=16)
    reference_y = tables.Float32Col(pos=17)
    reference_zenith = tables.Float32Col(pos=18)
    reference_azimuth = tables.Float32Col(pos=19)
    reference_size = tables.Float32Col(pos=20)
    reference_energy = tables.Float32Col(pos=21)


class ReconstructedEvent(ReconstructedCoincidence):

    """Store information about reconstructed events

    .. attribute:: id

        Index referring to the id of the event that was reconstructed.

    .. attribute:: d1,d2,d3,d4

        Booleans indicating which detectors participated in the
        reconstruction.

    """

    d1 = tables.BoolCol(pos=22)
    d2 = tables.BoolCol(pos=23)
    d3 = tables.BoolCol(pos=24)
    d4 = tables.BoolCol(pos=25)


class KascadeEvent(tables.IsDescription):

    """Store events from KASCADE"""

    run_id = tables.IntCol(pos=0)
    event_id = tables.Int64Col(pos=1)
    timestamp = tables.Time32Col(pos=2)
    nanoseconds = tables.UInt32Col(pos=3)
    ext_timestamp = tables.UInt64Col(pos=4)

    energy = tables.FloatCol(pos=5)
    core_pos = tables.FloatCol(pos=6, shape=2)
    zenith = tables.FloatCol(pos=7)
    azimuth = tables.FloatCol(pos=8)
    Num_e = tables.FloatCol(pos=9)
    Num_mu = tables.FloatCol(pos=10)
    dens_e = tables.FloatCol(pos=11, shape=4)
    dens_mu = tables.FloatCol(pos=12, shape=4)
    P200 = tables.FloatCol(pos=13)
    T200 = tables.FloatCol(pos=14)


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
