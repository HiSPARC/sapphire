""" Store CORSIKA simulation data in HDF5 file

    This module reads the CORSIKA binary ground particles file and stores
    each particle individually in a HDF5 file, using PyTables.  This file
    can then be used as input for the detector simulation.

"""
import argparse

import tables
from progressbar import ProgressBar, ETA, Bar, Percentage

from sapphire.corsika.reader import CorsikaFile


class GroundParticles(tables.IsDescription):
    """Store information about shower particles reaching ground level"""

    particle_id = tables.UInt16Col(pos=0)
    r = tables.Float32Col(pos=1)
    phi = tables.Float32Col(pos=2)
    x = tables.Float32Col(pos=3)
    y = tables.Float32Col(pos=4)
    t = tables.Float32Col(pos=5)
    p_x = tables.Float32Col(pos=6)
    p_y = tables.Float32Col(pos=7)
    p_z = tables.Float32Col(pos=8)
    hadron_generation = tables.UInt8Col(pos=9)
    observation_level = tables.UInt8Col(pos=10)


def save_particle(row, p):
    """Write the information of a particle into a row"""

    (p_x, p_y, p_z, x, y, t, id, r, hadron_generation, observation_level,
     phi) = p

    row['particle_id'] = id
    row['r'] = r
    row['phi'] = phi
    row['x'] = x
    row['y'] = y
    row['t'] = t
    row['p_x'] = p_x
    row['p_y'] = p_y
    row['p_z'] = p_z
    row['hadron_generation'] = hadron_generation
    row['observation_level'] = observation_level
    row.append()


def store_corsika_data(source, destination, table_name='groundparticles',
                       progress=False):
    """Store particles from a CORSIKA simulation in a HDF5 file

    :param source: CorsikaFile instance of the source DAT file
    :param destination: PyTables file instance of the destination file

    """
    if progress:
        print "Storing CORSIKA data (%s) in %s" % (source._filename,
                                                   destination.filename)
    source.check()

    for event in source.get_events():
        n_particles = event.get_end().n_particles_levels
        progress = progress and n_particles > 1
        try:
            table = destination.create_table('/', table_name, GroundParticles,
                                             'All groundparticles',
                                             expectedrows=n_particles)
        except tables.NodeError:
            if progress:
                print '%s already exists, doing nothing' % table_name
                return
            else:
                raise
        if progress:
            pbar = ProgressBar(maxval=n_particles - 1,
                               widgets=[Percentage(), Bar(), ETA()]).start()

        particle_row = table.row
        for row, particle in enumerate(event.get_particles()):
            save_particle(particle_row, particle)
            if progress and not row % 5000:
                pbar.update(row)
            if not row % 1000000:
                table.flush()

        if progress:
            pbar.finish()

    table.flush()

    run_header = source.get_header()
    run_end = source.get_end()
    for event in source.get_events():
        event_header = event.get_header()
        event_end = event.get_end()

    destination.set_node_attr('/', 'run_header', run_header)
    destination.set_node_attr('/', 'run_end', run_end)
    destination.set_node_attr('/', 'event_header', event_header)
    destination.set_node_attr('/', 'event_end', event_end)


def create_index(hdf_data, table_name='groundparticles', progress=False):
    """Create a completely sorted index for the x column

    This can speed up queries to select data based on the x column.

    """
    table = hdf_data.get_node('/', table_name)
    if progress:
        print 'Ensuring the x column for table %s is indexed.' % table_name
    try:
        table.cols.x.create_csindex()
    except ValueError:
        table.reindex_dirty()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', help="path of the CORSIKA source file")
    parser.add_argument('destination',
                        help="path of the HDF5 destination file")
    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite destination file it is already exists')
    parser.add_argument('--progress', action='store_true',
                        help='show progressbar during conversion')
    args = parser.parse_args()

    corsika_data = CorsikaFile(args.source)
    if args.overwrite:
        mode = 'w'
    else:
        mode = 'a'
    with tables.open_file(args.destination, mode) as hdf_data:
        store_corsika_data(corsika_data, hdf_data, progress=args.progress)
    with tables.open_file(args.destination, 'a') as hdf_data:
        create_index(hdf_data, progress=args.progress)


if __name__ == '__main__':
    main()
