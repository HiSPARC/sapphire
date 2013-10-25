#!/usr/bin/env python

import os
import sys
import shutil
import random
import textwrap
import subprocess

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
    THETAP    0.   0.                  range of zenith angle (degree)
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

    :param energy: the energy of the primary particle in log10(GeV)
    :param particle: string name of primary particle as given in particles
    :param queue: choose a queue to sumbit the job to:
                  short - max 4 hours, ...
                  stbcq - max 8 hours, 1000 jobs, 240 cpus
                  qlong - max 48 hours, 80 jobs
    :param corsika: name of the compiled CORSIKA executable to use:
                    corsika74000Linux_QGSII_gheisha
                    corsika74000Linux_EPOS_gheisha

    """
    def __init__(self, energy=7, particle='proton', queue='stbcq',
                 corsika='corsika74000Linux_QGSII_gheisha'):
        self.energy = energy
        self.particle = eval('particles.' + particle)
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

        inputpath = TEMPDIR + self.rundir + 'input-hisparc'
        input = INPUT_TEMPLATE.format(seed1=self.seed1, seed2=self.seed2,
                                      particle=self.particle,
                                      energy=self.energy)
        file = open(inputpath, 'w')
        file.write(input)
        file.close()

    def create_script(self):
        """Make Stoomboot script file"""

        scriptpath = TEMPDIR + self.rundir + 'run.sh'
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
        subprocess.check_output('ln -s {source}/* {dest}'
                                .format(source=CORSIKADIR,
                                        dest=TEMPDIR + self.rundir),
                                shell=True)


def multiple_jobs(n, E, p, q):
    """Use this to sumbit multiple jobs to stoomboot

    :param n: Number of jobs to submit
    :param E: log10(GeV) energy of primary particle
    :param p: particle kind (as string, see particles for possibilities)
    :param q: Stoomboot queue to submit to

    """
    for _ in arange(n):
        batch = CorsikaBatch(E, p, q)
        batch.batch_run()


if __name__ == '__main__':
    multiple_jobs(100, 7, 'proton', 'stbcq')
