""" Access the Nikhef Stoomboot cluster.

    .. note::
        This module is only for use at Nikhef. The Stoomboot cluster is only
        accessible for Nikhef users.

    Easy to use functions to make use of the Nikhef Stoomboot facilities.
    This checks the available slots on the requested queue, creates the
    scripts to submit, submits the jobs, and cleans up afterwards.

    Example usage::

        >>> from sapphire import qsub
        >>> qsub.check_queue('long')
        340
        >>> qsub.submit_job('touch /data/hisparc/test', 'job_1', 'express')

"""
import os
import subprocess

from . import utils


def check_queue(queue):
    """Check for available job slots on the selected queue for current user

    Maximum numbers from:
    ``qstat -Q -f | grep -e Queue: -e max_user_queuable -e max_queuable``.
    Note that some queues also have global maximum number of jobs.

    :param queue: queue name for which to check current number of job
                  slots in use.
    :return: number of available slots.

    """
    utils.which('qstat')
    all_jobs = int(subprocess.check_output('qstat {queue} | '
                                           'grep " [QR] " | wc -l'
                                           .format(queue=queue), shell=True))
    user_jobs = int(subprocess.check_output('qstat -u $USER {queue} | '
                                            'grep " [QR] " | wc -l'
                                            .format(queue=queue), shell=True))

    if queue == 'express':
        return 2 - user_jobs
    elif queue == 'short':
        return 1000 - user_jobs
    elif queue == 'generic':
        return min(2000 - user_jobs, 4000 - all_jobs)
    elif queue == 'long':
        return min(500 - user_jobs, 1000 - all_jobs)
    else:
        raise KeyError('Unknown queue name: {queue}'.format(queue=queue))


def submit_job(script, name, queue, extra=''):
    """Submit a job to Stoomboot

    :param script: contents of the script to run.
    :param name: name for the job.
    :param queue: name of the queue to run the job on.
    :param extra: optional extra arguments for the qsub command.

    """
    utils.which('qsub')
    script_path, script_name = create_script(script, name)

    # Effect of the arguments for qsub:
    # -q: the name of the queue to which the job is submitted
    # -V: all of the environment variables of the process are exported to
    #     the context of the batch job (e.g. PATH)
    # -z: do not print the job_identifier of the created job
    # -j oe: merge standard error into the standard output
    # -N: a recognizable name for the job
    qsub = ('qsub -q {queue} -V -z -j oe -N {name} {extra} {script}'
            .format(queue=queue, name=script_name, script=script_path,
                    extra=extra))

    result = subprocess.check_output(qsub, stderr=subprocess.STDOUT,
                                     shell=True)
    if not result == '':
        raise Exception('%s - Error occured: %s' % (name, result))

    delete_script(script_path)


def create_script(script, name):
    """Create script as temp file to be run on Stoomboot"""

    script_name = 'his_{name}.sh'.format(name=name)
    script_path = os.path.join('/tmp', script_name)

    with open(script_path, 'w') as script_file:
        script_file.write(script)
    os.chmod(script_path, 0774)

    return script_path, script_name


def delete_script(script_path):
    """Delete script after submitting to Stoomboot

    :param script_path: path to the script to be removed

    """
    os.remove(script_path)
