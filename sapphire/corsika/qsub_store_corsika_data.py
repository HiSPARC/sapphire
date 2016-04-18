""" Convert CORSIKA stored showers to HDF5 on Stoomboot

    Automatically submits Stoomboot jobs to convert corsika data. The
    script ``store_corsika_data`` can be used to convert a DAT000000
    CORSIKA file to a HDF5 file. This script checks our data folder for
    new or unconverted simulations and creates Stoomboot jobs to perform
    the conversion.

    This job is run as a cron job to ensure the simulations remain up to
    date.

"""
import os
import glob
import textwrap
import logging
import argparse

from .. import qsub


BIN_PATH = '/data/hisparc/env/miniconda/envs/corsika/bin/'
LOGFILE = '/data/hisparc/corsika/logs/qsub_store_corsika.log'
DATADIR = '/data/hisparc/corsika/data'
QUEUED_SEEDS = '/data/hisparc/corsika/queued.log'
SOURCE_FILE = 'DAT000000'
DESTINATION_FILE = 'corsika.h5'
SCRIPT_TEMPLATE = textwrap.dedent("""\
    #!/usr/bin/env bash
    umask 002
    {command}
    touch {datadir}
    # To alleviate Stoomboot, make sure the job is not to short.
    sleep $[ ( $RANDOM % 60 ) + 60 ]""")

logger = logging.getLogger('qsub_store_corsika_data')


def all_seeds():
    """Get set of all seeds in the corsika data directory"""

    dirs = glob.glob(os.path.join(DATADIR, '*_*'))
    seeds = [os.path.basename(dir) for dir in dirs]
    return set(seeds)


def seeds_processed():
    """Get the seeds of simulations for which the h5 is already created"""

    files = glob.glob(os.path.join(DATADIR, '*_*/corsika.h5'))
    seeds = [os.path.basename(os.path.dirname(file)) for file in files]
    return set(seeds)


def seeds_in_queue():
    """Get set of seeds already queued to be processed"""

    try:
        with open(QUEUED_SEEDS, 'r') as queued_seeds:
            seeds = queued_seeds.read().split('\n')
    except IOError:
        seeds = []
    return set(seeds)


def write_queued_seeds(seeds):
    """Write queued seeds to file"""

    with open(QUEUED_SEEDS, 'w') as queued_seeds:
        queued_seeds.write('\n'.join(seeds))


def append_queued_seeds(seeds):
    """Add seed to the queued seeds file"""

    queued = seeds_in_queue()
    updated_queued = queued.union(seeds).difference([''])
    write_queued_seeds(updated_queued)


def get_seeds_todo():
    """Get seeds to be processed"""

    seeds = all_seeds()
    processed = seeds_processed()
    queued = seeds_in_queue()

    queued = queued.difference(processed).difference([''])
    write_queued_seeds(queued)

    return seeds.difference(processed).difference(queued)


def filter_large_seeds(seeds_todo):
    """Exclude seeds for data files that are to large"""

    limit = 70e9  # larger than 70 GB has not been tested yet
    return {s for s in seeds_todo
            if os.path.getsize(os.path.join(DATADIR, s, SOURCE_FILE)) < limit}


def store_command(seed):
    """Write queued seeds to file"""

    source = os.path.join(DATADIR, seed, SOURCE_FILE)
    destination = os.path.join(DATADIR, seed, DESTINATION_FILE)
    command = ('{bin_path}python {bin_path}store_corsika_data {source} '
               '{destination}'.format(bin_path=BIN_PATH, source=source,
                                      destination=destination))

    return command


def run(queue):
    """Get list of seeds to process, then submit jobs to process them"""

    os.umask(002)
    logger.info('Getting todo list of seeds to convert.')
    seeds = get_seeds_todo()
    # seeds = filter_large_seeds(seeds)
    n_jobs_to_submit = min(len(seeds), qsub.check_queue(queue), 50)
    extra = ''
    if queue == 'long':
        extra += " -l walltime=96:00:00"

    logger.info('Submitting jobs for %d simulations.' % n_jobs_to_submit)
    try:
        for _ in xrange(n_jobs_to_submit):
            seed = seeds.pop()
            command = store_command(seed)
            script = SCRIPT_TEMPLATE.format(command=command, datadir=DATADIR)
            logger.info('Submitting job for %s.' % seed)
            qsub.submit_job(script, seed, queue, extra)
            append_queued_seeds([seed])
    except KeyError:
        logger.error('Out of seeds!')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-q', '--queue', metavar='name',
                        help="name of the Stoomboot queue to use, choose from "
                             "express, short, generic, and long (default)",
                        default='long',
                        choices=['express', 'short', 'generic', 'long'])
    args = parser.parse_args()
    logger.debug('Starting to submit new jobs.')
    run(args.queue)
    logger.info('Finished submitting jobs.')


if __name__ == '__main__':
    logging.basicConfig(
        filename=LOGFILE, filemode='a',
        format='%(asctime)s %(name)s %(levelname)s: %(message)s',
        datefmt='%y%m%d_%H%M%S', level=logging.INFO)
    main()
