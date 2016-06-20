""" Convert unsorted CORSIKA HDF5 to sorted HDF5 using Stoomboot"""

import os
import glob
import textwrap
import subprocess
import logging


QUEUE = 'generic'
LOGFILE = '/data/hisparc/corsika/logs/qsub_sort_corsika.log'
DATADIR = '/data/hisparc/corsika/data'
QUEUED_SEEDS = '/data/hisparc/corsika/sort_queued.log'
SCRIPT_TEMPLATE = textwrap.dedent("""\
    #!/usr/bin/env bash
    umask 002
    source activate corsika &> /dev/null
    cd {seed}
    ptrepack --sortby x --propindexes corsika.h5 corsika_sorted.h5
    mv corsika_sorted.h5 corsika.h5
    touch tmp_sorted_flag
    # To alleviate Stoomboot, make sure the job is not to short.
    sleep $[ ( $RANDOM % 60 ) + 60 ]""")

logging.basicConfig(filename=LOGFILE, filemode='a',
                    format='%(asctime)s %(name)s %(levelname)s: %(message)s',
                    datefmt='%y%m%d_%H%M%S', level=logging.INFO)
logger = logging.getLogger('qsub_store_corsika_data')


def all_seeds():
    """Get set of all seeds in the corsika data directory"""

    files = glob.glob(os.path.join(DATADIR, '*_*/corsika.h5'))
    seeds = [os.path.basename(os.path.dirname(file)) for file in files]
    return set(seeds)


def seeds_processed():
    """Get the seeds of simulations for which the h5 is already sorted"""

    files = glob.glob(os.path.join(DATADIR, '*_*/tmp_sorted_flag'))
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

    # Remove processed seeds and empty lines from queued log.
    queued = queued.difference(processed).difference([''])
    write_queued_seeds(queued)

    # Get seeds not yet processed and not in queue.
    return seeds.difference(processed).difference(queued)



def get_script_path(seed):
    """Create path for script"""

    script_name = 'sort_{seed}.sh'.format(seed=seed)
    script_path = os.path.join('/tmp', script_name)
    return script_path


def create_script(seed):
    """Create script as temp file to run on Stoomboot"""

    script_path = get_script_path(seed)
    input = SCRIPT_TEMPLATE.format(seed=os.path.join(DATADIR, seed))

    with open(script_path, 'w') as script:
        script.write(input)
    os.chmod(script_path, 0774)

    return script_path


def delete_script(seed):
    """Delete script"""

    script_path = get_script_path(seed)
    os.remove(script_path)


def submit_job(seed):
    """Submit job to Stoomboot"""

    script_path = create_script(seed)

    qsub = ('qsub -q {queue} -V -z -j oe -N {name} {script}'
            .format(queue=QUEUE, name=os.path.basename(script_path),
                    script=script_path))

    result = subprocess.check_output(qsub, stderr=subprocess.STDOUT,
                                     shell=True)
    if not result == '':
        msg = '%s - Error occured: %s' % (seed, result)
        logger.error(msg)
        raise Exception(msg)

    delete_script(seed)


def check_queue():
    """Check for available job slots on the chosen queue

    Maximum numbers from ``qstat -Q -f``

    :return: Number of available slots in the queue.

    """
    queued = 'qstat {queue} | grep [RQ] | wc -l'.format(queue=QUEUE)
    user_queued = ('qstat -u $USER {queue} | grep [RQ] | wc -l'
                   .format(queue=QUEUE))
    n_queued = int(subprocess.check_output(queued, shell=True))
    n_queued_user = int(subprocess.check_output(user_queued, shell=True))
    max_queue = 4000
    max_queue_user = 2000
    keep_free = 50

    return min(max_queue - n_queued,
               max_queue_user - n_queued_user) - keep_free


def run():
    """Get list of seeds to process, then submit jobs to process them"""

    os.umask(002)
    logger.info('Getting todo list of seeds to convert.')
    seeds = get_seeds_todo()
    n_jobs_to_submit = min(len(seeds), check_queue())
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
    logger.debug('Starting to submit new jobs.')
    run()
    logger.info('Finished submitting jobs.')


if __name__ == '__main__':
    main()
