""" Access the Nikhef Stoomboot cluster.

    .. note::
        This module is only for use at Nikhef. The Stoomboot cluster is only
        accessible for Nikhef users.

    Easy to use functions to make use of the Nikhef Stoomboot facilities.
    This checks the available slots on the requested queue, creates the
    scripts to submit, submits the jobs, and cleans up afterwards.

    Example usage::

        >>> from sapphire import nsub
        >>> nsub.submit_job('touch /data/hisparc/test', 'job_1', 'express')

"""
import os
import subprocess

def submit_job(script, name, queue, extra=''):
    """Submit a job to Stoomboot

    :param script: contents of the script to run.
    :param name: name for the job.
    :param queue: name of the queue to run the job on.
    :param extra: optional extra arguments for the nsub command.

    """
    script_path, script_name = create_script(script, name)

    # Effect of the arguments for nsub:
    # -q: the name of the queue to which the job is submitted
    # -V: all of the environment variables of the process are exported to
    #     the context of the batch job (e.g. PATH)
    # -z: do not print the job_identifier of the created job
    # -j oe: merge standard error into the standard output
    # -N: a recognizable name for the job
    nsub = ('/global/ices/toolset/bin/nsub -q {queue} -N {name} {extra} {script}'
            .format(queue=queue, name=script_name, script=script_path, extra=extra))

    result = subprocess.check_output(nsub, stderr=subprocess.STDOUT,
                                     shell=True)
    if not result == b'':
        raise Exception('%s - Error occured: %s' % (name, result))

    delete_script(script_path)


def create_script(script, name):
    """Create script as temp file to be run on Stoomboot"""

    script_name = 'his_{name}.sh'.format(name=name)
    script_path = os.path.join('/data/tunnel/user/kaspervd/tmp', script_name)

    with open(script_path, 'w') as script_file:
        script_file.write(script)
    os.chmod(script_path, 0o774)

    return script_path, script_name


def delete_script(script_path):
    """Delete script after submitting to Stoomboot

    :param script_path: path to the script to be removed

    """
    os.remove(script_path)
