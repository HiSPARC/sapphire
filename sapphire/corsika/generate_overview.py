import os
import os.path
import tables
import glob
import logging

import progressbar as pb

from sapphire import corsika

LOGFILE = '/data/hisparc/corsika/logs/generate_overview.log'
DAT_URL = '/data/hisparc/corsika/data'
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


def save_seed(row, seeds, header, end):
    """Write the information of one simulation into a row

    :param row: a new row instance to be appended to the table
    :param seeds: the unique id consisting of the two seeds
    :param header, end: the event header and end for the simulation

    """
    seed1, seed2 = seeds.split('_')
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


def write_row(output_row, seeds):
    """Read the header of a simulation and write this to the output."""

    try:
        with tables.openFile(os.path.join(DAT_URL, seeds, 'corsika.h5'), 'r') \
                as corsika_data:
            try:
                groundparticles = corsika_data.getNode('/groundparticles')
                header = groundparticles._v_attrs.event_header
                end = groundparticles._v_attrs.event_end
                save_seed(output_row, seeds, header, end)
            except tables.NoSuchNodeError:
                logger.INFO('No groundparticles table for %s' % seeds)
            except AttributeError:
                logger.INFO('Missing attribute (header or footer) for %s' %
                            seeds)
    except tables.HDF5ExtError:
        logger.INFO('Unable to open file for %s' % seeds)


def get_simulations(simulations_data):
    """Get the information of the simulations and create a table."""

    files = glob.glob(os.path.join(DAT_URL, '*/corsika.h5'))
    simulations_table = simulations_data.getNode('/simulations')
    progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])
    for file in progress(files):
        output_row = simulations_table.row
        dir = os.path.dirname(file)
        seeds = os.path.basename(dir)
        write_row(output_row, seeds)
        simulations_table.flush()


def prepare_output():
    """Write the table to seed_info.h5"""

    simulations_data = tables.openFile(os.path.join(OUTPUT_PATH,
                                                    'seed_info.h5'), 'w')
    simulations_data.createTable('/', 'simulations', Simulations,
                                 'Simulations overview')
    return simulations_data


def generate_simulation_overview():
    simulations_data = prepare_output()
    get_simulations(simulations_data)
    simulations_data.close()


if __name__ == '__main__':
    generate_simulation_overview()
