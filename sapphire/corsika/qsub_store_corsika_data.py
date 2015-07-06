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
import subprocess
import logging
import argparse


LOGFILE = '/data/hisparc/corsika/logs/qsub_store_corsika.log'
DATADIR = '/data/hisparc/corsika/data'
QUEUED_SEEDS = '/data/hisparc/corsika/queued.log'
SOURCE_FILE = 'DAT000000'
DESTINATION_FILE = 'corsika.h5'
SCRIPT_TEMPLATE = textwrap.dedent("""\
    #!/usr/bin/env bash
    umask 002
    source activate /data/hisparc/corsika_env &> /dev/null
    {command}
    touch {datadir}
    # To alleviate Stoomboot, make sure the job is not to short.
    sleep $[ ( $RANDOM % 60 ) + 60 ]""")

logging.basicConfig(filename=LOGFILE, filemode='a',
                    format='%(asctime)s %(name)s %(levelname)s: %(message)s',
                    datefmt='%y%m%d_%H%M%S', level=logging.INFO)
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
        for seed in seeds:
            queued_seeds.write('%s\n' % seed)


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


def store_command(seed):
    """Write queued seeds to file"""

    source = os.path.join(DATADIR, seed, SOURCE_FILE)
    destination = os.path.join(DATADIR, seed, DESTINATION_FILE)
    command = 'store_corsika_data %s %s' % (source, destination)

    return command


def create_script(seed):
    """Create script as temp file to run on Stoomboot"""

    script_name = 'store_{seed}.sh'.format(seed=seed)
    script_path = os.path.join('/tmp', script_name)
    command = store_command(seed)
    input = SCRIPT_TEMPLATE.format(command=command, datadir=DATADIR)

    with open(script_path, 'w') as script:
        script.write(input)
    os.chmod(script_path, 0774)

    return script_path, script_name


def delete_script(seed):
    """Delete script"""

    script_name = 'store_{seed}.sh'.format(seed=seed)
    script_path = os.path.join('/tmp', script_name)

    os.remove(script_path)


def submit_job(seed):
    """Submit job to Stoomboot"""

    script_path, script_name = create_script(seed)

    qsub = ('qsub -q {queue} -V -z -j oe -N {name} {script}'
            .format(queue=QUEUE, name=script_name, script=script_path))

    result = subprocess.check_output(qsub, stderr=subprocess.STDOUT,
                                     shell=True)
    if not result == '':
        logger.error('%s - Error occured: %s' % (seed, result))
        raise Exception

    delete_script(seed)


def check_queue(queue):
    """Check for available job slots on the selected queue for current user

    Maximum numbers from ``qstat -Q -f``

    :param queue: queue name for which to check current number of job
                  slots in use.
    :return: boolean, True if slots are available, False otherwise.

    """
    all_jobs = int(subprocess.check_output('qstat {queue} | '
                                           'grep " [QR] " | wc -l'
                                           .format(queue=QUEUE), shell=True))
    user_jobs = int(subprocess.check_output('qstat -u $USER {queue} | '
                                            'grep " [QR] " | wc -l'
                                            .format(queue=QUEUE), shell=True))

    if QUEUE == 'express':
        return 2 - user_jobs
    elif QUEUE == 'short':
        return 1000 - user_jobs
    elif QUEUE == 'generic':
        return min(2000 - user_jobs, 4000 - all_jobs)
    elif QUEUE == 'long':
        return min(500 - user_jobs, 1000 - all_jobs)
    else:
        raise KeyError('Unknown queue name: {queue}'.format(queue=QUEUE))


def run(queue):
    """Get list of seeds to process, then submit jobs to process them"""

    os.umask(002)
    logger.info('Getting todo list of seeds to convert.')
    seeds = get_seeds_todo()
    n_jobs_to_submit = min(len(seeds), check_queue(queue))
    logger.info('Submitting jobs for %d simulations.' % n_jobs_to_submit)
    try:
        for _ in xrange(n_jobs_to_submit):
            seed = seeds.pop()
            logger.info('Submitting job for %s.' % seed)
            submit_job(seed)
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
    main()
