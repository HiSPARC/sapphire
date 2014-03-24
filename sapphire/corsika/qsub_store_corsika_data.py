import os
import glob
import textwrap
import subprocess


DATADIR = '/data/hisparc/corsika/data'
QUEUED_SEEDS = '/data/hisparc/corsika/queued.log'
SOURCE_FILE = 'DAT000000'
DESTINATION_FILE = 'corsika.h5'
SCRIPT_TEMPLATE = textwrap.dedent("""\
    #!/usr/bin/env bash
    umask 002
    {command}""")


def all_seeds():
    """Get set of all seeds in the corsika data directory"""

    dirs = glob.glob(os.path.join(DATADIR, '*_*'))
    seeds = [os.path.basename(dir) for dir in dirs]
    return set(seeds)


def seeds_processed():
    """Get the seeds for the simualtions for with an h5 is already created"""

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
    updated_queued = queued.union(seeds)
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
    input = SCRIPT_TEMPLATE.format(command=command)

    with open(script_path, 'w') as scipt:
        scipt.write(input)
    os.chmod(script_path, 0774)

    return script_path, script_name

def submit_job(seed):
    """Submit job to Stoomboot"""

    script_path, script_name = create_script(seed)

    qsub = ('qsub -q short -V -z -j oe -N {name} {script}'
            .format(name=script_name, script=script_path))

    try:
        result = subprocess.check_output(qsub, stderr=subprocess.STDOUT,
                                         shell=True)
    if not result == '':
        print '%s - Error occured: %s' % (self.seed, result)
        raise Exception


def run():
    """Get list of seeds to process, then submit jobs to process them"""

    seeds = get_seeds_todo()
    for seed in seeds[:10]:
        submit_job(seed)
        append_queued_seeds([seed])


if __name__ == '__main__':
    run()
