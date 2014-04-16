import os
import tables
import glob
import logging
import shutil
import argparse

import progressbar as pb

from sapphire import corsika

LOGFILE = '/data/hisparc/corsika/logs/generate_overview.log'
DATA_PATH = '/data/hisparc/corsika/data'
OUTPUT_PATH = '/data/hisparc/corsika'

logging.basicConfig(filename=LOGFILE, filemode='a',
                    format='%(asctime)s %(name)s %(levelname)s: %(message)s',
                    datefmt='%y%m%d_%H%M%S', level=logging.INFO)
logger = logging.getLogger('generate_overview')


class Simulations(tables.IsDescription):
    """Store information about shower particles reaching ground level"""

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
    :param header, end: the event header and end for the simulation.

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


def read_seeds(simulations_table, seeds):
    """Read the header of a simulation and write this to the output."""

    try:
        with tables.openFile(os.path.join(DATA_PATH, seeds, 'corsika.h5'),
                             'r') as corsika_data:
            try:
                header = corsika_data.root._v_attrs.event_header
                end = corsika_data.root._v_attrs.event_end
                write_row(simulations_table, seeds, header, end)
            except AttributeError:
                logger.info('%19s: Missing attribute (header or end).' % seeds)
    except (IOError, tables.HDF5ExtError):
        logger.info('%19s: Unable to open file.' % seeds)


def get_simulations(simulations, overview, progress=False):
    """Get the information of the simulations and create a table."""

    simulations_table = overview.getNode('/simulations')
    if progress:
        pbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])
        simulations = pbar(simulations)
    for seeds in simulations:
        read_seeds(simulations_table, seeds)
    simulations_table.flush()


def prepare_output(n):
    """Create a temporary file in which to store the overview

    :param n: the number of simulations, i.e. expected number of rows.

    """
    os.umask(002)
    tmp_output = os.path.join(OUTPUT_PATH, 'temp_simulation_overview.h5')
    overview = tables.openFile(tmp_output, 'w')
    overview.createTable('/', 'simulations', Simulations,
                         'Simulations overview', expectedrows=n)
    return overview


def move_tempfile_to_destination():
    tmp_path = os.path.join(OUTPUT_PATH, 'temp_simulation_overview.h5')
    data_path = os.path.join(OUTPUT_PATH, 'simulation_overview.h5')
    shutil.move(tmp_path, data_path)


def generate_simulation_overview(progress=False):
    logger.info('Getting simulation list.')
    simulations = os.walk(DATA_PATH).next()[1]
    overview = prepare_output(len(simulations))
    get_simulations(simulations, overview, progress=progress)
    overview.close()
    move_tempfile_to_destination()
    logger.info('Finished generating overview.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--progress', action='store_true',
                        help='show progressbar during generation')
    args = parser.parse_args()
    generate_simulation_overview(progress=args.progress)
