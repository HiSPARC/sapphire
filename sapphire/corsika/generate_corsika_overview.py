""" Generate and overview table of the CORSIKA simulations

    This script will look for all completed and converted CORSIKA
    simulations in the given data path. Information about each
    simulation is collected and then summarized in a new h5 file as an
    overview.

    The given source path should contain subdirectories named after the
    seeds used for the simulation in the format ``{seed1}_{seed2}``,
    e.g. ``821280921_182096636``. These in turn should contain converted
    CORSIKA simulation results called ``corsika.h5``.

"""
import os
import glob
import tables
import logging
import shutil
import argparse
import tempfile

from ..utils import pbar


LOGFILE = '/data/hisparc/corsika/logs/generate_overview.log'
DATA_PATH = '/data/hisparc/corsika/data'
OUTPUT_PATH = '/data/hisparc/corsika/corsika_overview.h5'

logger = logging.getLogger('generate_corsika_overview')


class Simulations(tables.IsDescription):
    """Store summary information about CORSIKA simulations"""

    seed1 = tables.UInt32Col(pos=0)
    seed2 = tables.UInt32Col(pos=1)
    particle_id = tables.UInt32Col(pos=2)
    energy = tables.Float32Col(pos=3)
    first_interaction_altitude = tables.Float32Col(pos=4)
    p_x = tables.Float32Col(pos=5)
    p_y = tables.Float32Col(pos=6)
    p_z = tables.Float32Col(pos=7)
    zenith = tables.Float32Col(pos=8)
    azimuth = tables.Float32Col(pos=9)
    observation_height = tables.Float32Col(pos=10)
    n_photon = tables.Float32Col(pos=11)
    n_electron = tables.Float32Col(pos=12)
    n_muon = tables.Float32Col(pos=13)
    n_hadron = tables.Float32Col(pos=14)


def write_row(table, seeds, header, end):
    """Write the information of one simulation into a row

    :param table: the table where the new data should be appended.
    :param seeds: the unique id consisting of the two seeds.
    :param header,end: the event header and end for the simulation.

    """
    seed1, seed2 = seeds.split('_')
    row = table.row
    row['seed1'] = seed1
    row['seed2'] = seed2
    row['particle_id'] = header.particle_id
    row['energy'] = header.energy
    row['first_interaction_altitude'] = header.first_interaction_altitude
    row['p_x'] = header.p_x
    row['p_y'] = header.p_y
    row['p_z'] = header.p_z
    row['zenith'] = header.zenith
    row['azimuth'] = header.azimuth
    row['observation_height'] = header.observation_heights[0]
    row['n_photon'] = end.n_photons_levels
    row['n_electron'] = end.n_electrons_levels
    row['n_muon'] = end.n_muons_levels
    row['n_hadron'] = end.n_hadrons_levels
    row.append()


def read_seeds(simulations_table, source, seeds):
    """Read the header and end of a simulation and write to overview.

    :param simulations_table: PyTables table in which the simulation
                              overview is stored.
    :param source: directory containing the CORSIKA simulations.
    :param seeds: directory name of a simulation, format: '{seed1}_{seed2}'.

    """
    path = os.path.join(source, seeds, 'corsika.h5')
    if not os.path.exists(path):
        logger.info('%19s: No corsika.h5 available.' % seeds)
        return
    try:
        with tables.open_file(path, 'r') as corsika_data:
            try:
                header = corsika_data.get_node_attr('/', 'event_header')
                end = corsika_data.get_node_attr('/', 'event_end')
                write_row(simulations_table, seeds, header, end)
            except AttributeError:
                logger.info('%19s: Missing attribute (header or end).' % seeds)
    except (IOError, tables.HDF5ExtError):
        logger.info('%19s: Unable to open file.' % seeds)


def get_simulations(source, simulations, overview, progress=False):
    """Get the information of the simulations and create a table."""

    simulations_table = overview.get_node('/simulations')
    for seeds in pbar(simulations, show=progress):
        read_seeds(simulations_table, source, seeds)
    simulations_table.flush()


def prepare_output(n):
    """Create a temporary file in which to store the overview

    :param n: the number of simulations, i.e. expected number of rows.
    :return: path to the temporary file and a PyTables handler for the file.

    """
    os.umask(002)
    tmp_path = create_tempfile_path()
    overview = tables.open_file(tmp_path, 'w')
    overview.create_table('/', 'simulations', Simulations,
                          'Simulations overview', expectedrows=n)
    return tmp_path, overview


def create_tempfile_path():
    fd, path = tempfile.mkstemp('.h5')
    os.close(fd)
    return path


def move_tempfile_to_destination(tmp_path, destination):
    shutil.move(tmp_path, destination)


def all_seeds(source):
    """Get set of all seeds in the corsika data directory"""

    dirs = glob.glob(os.path.join(source, '*_*'))
    seeds = [os.path.basename(dir) for dir in dirs]
    return set(seeds)


def generate_corsika_overview(source, destination, progress=False):
    logger.info('Getting simulation list.')
    # Get names of all subdirectories
    simulations = all_seeds(source)
    tmp_path, overview = prepare_output(len(simulations))
    get_simulations(source, simulations, overview, progress=progress)
    overview.close()
    move_tempfile_to_destination(tmp_path, destination)
    logger.info('Finished generating overview.')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('source', nargs='?', default=DATA_PATH,
                        help="directory path containing CORSIKA simulations")
    parser.add_argument('destination', nargs='?', default=OUTPUT_PATH,
                        help="path of the HDF5 output file")
    parser.add_argument('--progress', action='store_true',
                        help='show progressbar during generation')
    parser.add_argument('--log', action='store_true',
                        help='write logs to file, only for use on server')
    parser.add_argument('--lazy', action='store_true',
                        help='only run if the overview is outdated')
    args = parser.parse_args()
    if args.log:
        logging.basicConfig(filename=LOGFILE, filemode='a',
                            format='%(asctime)s %(name)s %(levelname)s: '
                                   '%(message)s',
                            datefmt='%y%m%d_%H%M%S', level=logging.INFO)
    if args.lazy:
        last_store = os.path.getmtime(args.source)
        last_overview = os.path.getmtime(args.destination)
        if last_overview > last_store:
            logger.info('Overview up to date.')
            return

    generate_corsika_overview(source=args.source,
                              destination=args.destination,
                              progress=args.progress)


if __name__ == '__main__':
    main()
