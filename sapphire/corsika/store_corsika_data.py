""" Store CORSIKA simulation data in HDF5 file

    This module reads the CORSIKA binary ground particles file and stores
    each particle individually in a HDF5 file, using PyTables.  This file
    can then be used as input for the detector simulation.

    The syntax and options for calling this script can be seen with::

        $ store_corsika_data --help

    For example to convert a CORSIKA file in the current directory called
    DAT000000 to a HDF5 called corsika.h5 with a progress bar run::

        $ store_corsika_data --progress DAT000000 corsika.h5

"""
import argparse
import tempfile
import os

import tables
from progressbar import ProgressBar, ETA, Bar, Percentage

from .reader import CorsikaFile
from .mergesort import TableMergeSort


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


def store_and_sort_corsika_data(source, destination, overwrite=False,
                                progress=False):
    """First convert the data to HDF5 and create a sorted version"""

    if os.path.exists(destination):
        if not overwrite:
            if progress:
                raise Exception("Destination already exists, doing nothing")
            return
        else:
            os.remove(destination)

    corsika_data = CorsikaFile(source)

    temp_dir = os.path.dirname(destination)
    unsorted = create_tempfile_path(temp_dir)
    temp_path = create_tempfile_path(temp_dir)

    with tables.open_file(unsorted, 'a') as hdf_temp:
        store_corsika_data(corsika_data, hdf_temp, progress=progress)
    with tables.open_file(unsorted, 'r') as hdf_unsorted, \
            tables.open_file(destination, 'w') as hdf_data, \
            tables.open_file(temp_path, 'w') as hdf_temp:

        with TableMergeSort('x', hdf_unsorted, hdf_data, hdf_temp,
                            progress=progress) as mergesort:
            mergesort.sort()

            event_header = hdf_unsorted.get_node_attr('/', 'event_header')
            run_header = hdf_unsorted.get_node_attr('/', 'run_header')
            event_end = hdf_unsorted.get_node_attr('/', 'event_end')
            run_end = hdf_unsorted.get_node_attr('/', 'run_end')
            hdf_data.set_node_attr('/', 'event_header', event_header)
            hdf_data.set_node_attr('/', 'run_header', run_header)
            hdf_data.set_node_attr('/', 'event_end', event_end)
            hdf_data.set_node_attr('/', 'run_end', run_end)

    os.remove(unsorted)
    os.remove(temp_path)

    with tables.open_file(destination, 'a') as hdf_data:
        create_index(hdf_data, progress=progress)


def store_corsika_data(source, destination, table_name='groundparticles',
                       progress=False):
    """Store particles from a CORSIKA simulation in a HDF5 file

    :param source: CorsikaFile instance of the source DAT file
    :param destination: PyTables file instance of the destination file

    """
    if progress:
        print "Converting CORSIKA data (%s) to HDF5 format" % source._filename
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


def copy_and_sort_node(hdf_temp, hdf_data, table_name='groundparticles',
                       progress=False):
    """Sort the data in the tables by the x column

    This speeds up queries to select data based on the x column.

    """
    target_root = hdf_data.get_node('/')
    source_table = hdf_temp.get_node('/', table_name)
    if progress:
        print 'Creating the sorted HDF5 file.'
    source_table.copy(newparent=target_root, sortby='x', propindexes=True)
    hdf_temp.copy_node_attrs('/', target_root)


def create_tempfile_path(temp_dir=None):
    """Create a temporary file, close it, and return the path"""

    f, path = tempfile.mkstemp(suffix='.h5', dir=temp_dir)
    os.close(f)
    return path


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

    store_and_sort_corsika_data(args.source, args.destination, args.overwrite,
                                args.progress)


if __name__ == '__main__':
    main()
