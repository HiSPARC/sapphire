#!/usr/bin/env python

import os
import sys
import shutil
import random
import textwrap
import subprocess

import progressbar as pb

import particles


TEMPDIR = '/data/hisparc/corsika/running/'
DATADIR = '/data/hisparc/corsika/data/'
CORSIKADIR = '/data/hisparc/corsika/corsika-74000/run/'
INPUT_TEMPLATE = textwrap.dedent("""\
    RUNNR     0                        run number
    EVTNR     1                        number of first shower event
    SEED      {seed1}   0   0          seed for 1. random number sequence (hadron shower)
    SEED      {seed2}   0   0          seed for 2. random number sequence (EGS4)
    SEED      3   0   0                seed for 3. random number sequence (Cherenkov)
    SEED      4   0   0                seed for 4. random number sequence (IACT)
    SEED      5   0   0                seed for 5. random number sequence (NUPRIM)
    SEED      6   0   0                seed for 6. random number sequence (PARALLEL)
    NSHOW     1                        number of showers to generate (MAX 1 for PARALLEL)
    PRMPAR    {particle}               particle type of prim. particle (14=proton, 1=photon, 3=electron)
    ERANGE    1.E{energy}  1.E{energy} energy range of primary particle (GeV)
    ESLOPE    -2.7                     slope of primary energy spectrum (E^y)
    THETAP    {theta}   {theta}        range of zenith angle (degree)
    PHIP      0.   0.                  range of azimuth angle (degree)
    FIXCHI    0.                       starting altitude (g/cm**2)
    FIXHEI    0.   0                   height and target type of first interaction (cm, [1=N, 2=O, 3=Ar])
    MAGNET    18.908 45.261            magnetic field DAin Amsterdam (uT)
    HADFLG    0  0  0  0  0  2         flags hadr.interact.&fragmentation
    ELMFLG    T   T                    em. interaction flags (NKG,EGS)
    STEPFC    1.0                      mult. scattering step length fact.
    RADNKG    1000.E2                  outer radius for NKG lat.dens.distr. (cm)
    ECUTS     0.3  0.3  0.003  0.003   energy cuts for particles (GeV, hadrons, muons, electrons, photons)
    LONGI     T  10.  T  T             longit.distr. & step size & fit & out
    MUMULT    T                        muon multiple scattering angle
    MUADDI    F                        additional info for muons
    OBSLEV    10.E2                    observation level -maaiveld -3.7m (cm)
    MAXPRT    1                        max. number of printed events
    ECTMAP    1.E4                     cut on gamma factor for printout
    DIRECT    ./                       output directory
    USER      hisparc                  user
    DEBUG     F  6  F  1000000         debug flag and log.unit for out
    EXIT                               terminates input""")
SCRIPT_TEMPLATE = textwrap.dedent("""\
    #!/usr/bin/env bash

    umask 002

    NAME="his_{seed1}_{seed2}"
    SCRIPT=/tmp/$NAME

    # Start Stoomboot script
    cat >> $SCRIPT << EOF
    #!/usr/bin/env bash

    umask 002

    echo "--- PBS job configuration ---"
    echo "PBS_O_HOST = " \$PBS_O_HOST
    echo "PBS_O_QUEUE = " \$PBS_O_QUEUE
    echo "PBS_O_WORKDIR = " \$PBS_O_WORKDIR
    echo "PBS_ENVIRONMENT = " \$PBS_ENVIRONMENT
    echo "PBS_JOBID = " \$PBS_JOBID
    echo "PBS_JOBNAME = " \$PBS_JOBNAME
    echo "PBS_NODEFILE = " \$PBS_NODEFILE
    echo "PBS_QUEUE = " \$PBS_QUEUE
    echo "-----------------------------"

    export PATH=\${{PBS_O_PATH}}

    # Run CORSIKA
    /usr/bin/time -o time.log ./{corsika} < input-hisparc > corsika-output.log

    # Clean up after run
    find -type l -delete
    cd ../..
    mv {rundir} data/

    EOF
    # End of Stoomboot script

    chmod ug+x $SCRIPT

    qsub -N ${{NAME}} -q {queue} -V -j oe -d {rundir} $SCRIPT

    rm $SCRIPT""")


class CorsikaBatch(object):
    """Run many simultaneous CORSIKA simulations using Stoomboot

    Stoomboot is the Nikhef computer cluster.

    :param energy: the energy of the primary particle in log10(GeV),
                   so an energy of 7 corresponds to 10**6 eV.
    :param particle: name of primary particle, e.g. proton, gamma or iron
    :param theta: zenith angle of the primary particle (in degrees),
                  common choices: 0, 7.5, 15, 22.5, 30, 37.5 and 45.
    :param queue: choose a queue to sumbit the job to:
                  short - max 4 hours, max 1000 jobs
                  stbcq - max 8 hours, max 1000 jobs
                  qlong - max 48 hours, max 500 jobs, max 64 running
    :param corsika: name of the compiled CORSIKA executable to use:
                    corsika74000Linux_QGSII_gheisha
                    corsika74000Linux_EPOS_gheisha

    """
    def __init__(self, energy=7, particle='proton', theta=22.5, queue='stbcq',
                 corsika='corsika74000Linux_QGSII_gheisha'):
        self.energy = energy
        self.particle = [i for i, p in particles.id.items()
                         if p == particle][0]
        self.theta = theta
        self.queue = queue
        self.corsika = corsika
        self.seed1 = None
        self.seed2 = None
        self.rundir = None

    def batch_run(self):
        self.prepare_env()
        self.submit_job()

    def prepare_env(self):
        """Setup CORSIKA environment"""

        # Set umask
        os.umask(002)

        # Setup directories
        taken = self.taken_seeds()
        self.gen_random_seeds(taken)
        self.make_rundir()
        self.goto_rundir()

        # Create/copy files
        self.create_input()
        self.create_script()
        self.symlink_corsika()

    def submit_job(self):
        """Submit job to Stoomboot"""

        subprocess.check_output('./run.sh', shell=True)

    def taken_seeds(self):
        """Get list of seeds already used"""

        taken = os.listdir(DATADIR)
        taken.extend(os.listdir(TEMPDIR))
        return taken

    def gen_random_seeds(self, taken):
        """Get unused combination of two seeds for CORSIKA

        :param taken: List of seed combinations already taken
                      each is formatted like this: 'seed1_seed2'

        """
        seed1 = random.randint(1, 900000000)
        seed2 = random.randint(1, 900000000)
        seed = "{seed1}_{seed2}".format(seed1=seed1, seed2=seed2)
        if not seed in taken:
            self.seed1 = seed1
            self.seed2 = seed2
            self.rundir = seed +'/'
        else:
            self.get_random_seeds(taken)

    def make_rundir(self):
        """Make the run directory"""

        os.mkdir(TEMPDIR + self.rundir)

    def goto_rundir(self):
        """Move into the run directory"""

        os.chdir(TEMPDIR + self.rundir)

    def create_input(self):
        """Make CORSIKA steering file"""

        inputpath = os.path.join(TEMPDIR, self.rundir, 'input-hisparc')
        input = INPUT_TEMPLATE.format(seed1=self.seed1, seed2=self.seed2,
                                      particle=self.particle,
                                      energy=self.energy, theta=self.theta)
        file = open(inputpath, 'w')
        file.write(input)
        file.close()

    def create_script(self):
        """Make Stoomboot script file"""

        scriptpath = os.path.join(TEMPDIR, self.rundir, 'run.sh')
        script = SCRIPT_TEMPLATE.format(seed1=self.seed1, seed2=self.seed2,
                                        queue=self.queue, corsika=self.corsika,
                                        rundir=TEMPDIR + self.rundir)
        file = open(scriptpath, 'w')
        file.write(script)
        file.close()
        os.chmod(scriptpath, 0774)

    def symlink_corsika(self):
        """Create symbolic links to CORSIKA run files

        CORSIKA requires files from the run directory to be available in
        the PWD. So we create symlinks to all files in the run dir.

        """
        source = os.path.join(CORSIKADIR, '*')
        destination = os.path.join(TEMPDIR, self.rundir)
        subprocess.check_output('ln -s {source} {destination}'
                                .format(source=source,
                                        destination=destination),
                                shell=True)


def check_queue(queue):
    """Check for available job slots on the selected queue for current user

    Caveat: also counts completed jobs still listed in qstat,
            add `grep [RQ]` to fix.

    :param queue: queue name for which to check current number of job
                  slots in use.
    :return: boolean, True if slots are available, False otherwise.

    """
    all_jobs = int(subprocess.check_output('qstat | grep {queue} | wc -l'
                                           .format(queue=queue), shell=True))
    user_jobs = int(subprocess.check_output('qstat -u $USER | grep {queue}'
                                            ' | wc -l'.format(queue=queue),
                                            shell=True))

    if queue == 'short':
        return 1000 - user_jobs
    elif queue == 'stbcq':
        return min(1000 - user_jobs, 10000 - all_jobs)
    elif queue == 'qlong':
        return min(500 - user_jobs, 1000 - all_jobs)
    else:
        raise KeyError('Unknown queue name: {queue}'.format(queue))


def multiple_jobs(n, E, p, T, q, c):
    """Use this to sumbit multiple jobs to Stoomboot

    :param n: Number of jobs to submit
    :param E: log10(GeV) energy of primary particle
    :param p: Particle kind (as string, see particles.py for possibilities)
    :param T: Zenith angle in degrees of the primary particle
    :param q: Stoomboot queue to submit to
    :param c: Name of the CORSIKA executable to use

    """
    print textwrap.dedent("""\
        Batch submitting jobs to Stoomboot:
        Number of jobs      {n}
        Particle energy     10^{E} GeV
        Primary particle    {p}
        Zenith angle        {T} degrees
        Stoomboot queue     {q}
        CORSIKA executable  {c}
        """.format(n=n, E=E, p=p, T=T, q=q, c=c))

    available_slots = check_queue(q)
    if available_slots <= 0:
        n = 0
        print 'Submitting no jobs because queue is full.'
        return
    elif available_slots < n:
        n = available_slots
        print 'Submitting {n} jobs because queue is almost full.'.format(n=n)

    progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(), pb.ETA()])
    for _ in progress(xrange(n)):
        batch = CorsikaBatch(energy=E, particle=p, theta=T, queue=q, corsika=c)
        batch.batch_run()
    print 'Done.'


if __name__ == '__main__':
    multiple_jobs(100, 7, 'proton', 22.5, 'stbcq',
                  'corsika74000Linux_QGSII_gheisha')
