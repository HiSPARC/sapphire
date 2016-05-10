""" Run CORSIKA simulations on Stoomboot

    In order to quickly get a good sample of simulated showers we use the
    Nikhef computer cluster Stoomboot to run multiple jobs simultaneously.
    For this purpose a script has been written that will make this easy.
    The :mod:`~sapphire.corsika.qsub_corsika` script can submit as many
    jobs as you want with the parameters that you desire. It automatically
    ensures that a unique combination of seeds for the random number
    sequences are used for each simulation.

    To run this file correctly do it in the correct env::

        $ source activate corsika

    The syntax for calling the script can be seen by calling its help::

        $ qsub_corsika --help

    For example, running 100 showers with proton primaries of 1e16 eV
    coming in at 22.5 degrees zenith and 90 degrees azimuth on the
    standard Stoomboot queue with the default CORSIKA configuration::

        $ qsub_corsika 100 16 proton 22.5 -q generic -a 90

"""
import os
import random
import textwrap
import subprocess
import argparse
import warnings
from math import modf

from . import particles
from ..utils import pbar
from .. import qsub


TEMPDIR = '/data/hisparc/corsika/running/'
DATADIR = '/data/hisparc/corsika/data/'
FAILDIR = '/data/hisparc/corsika/failed/'
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
    ERANGE    {energy_pre}E{energy_pow}  {energy_pre}E{energy_pow} energy range of primary particle (GeV)
    ESLOPE    -2.7                     slope of primary energy spectrum (E^y)
    THETAP    {theta}   {theta}        range of zenith angle (degree)
    PHIP      {phi}   {phi}            range of azimuth angle, phi is direction the shower points to (degree)
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
    DATDIR    {tablesdir}              location of the input data tables
    DIRECT    ./                       output directory
    USER      hisparc                  user
    DEBUG     F  6  F  1000000         debug flag and log.unit for out
    EXIT                               terminates input""")
SCRIPT_TEMPLATE = textwrap.dedent("""\
    #!/usr/bin/env bash

    umask 002

    # Run CORSIKA
    /usr/bin/time -o time.log {corsika} < input-hisparc > corsika-output.log

    # Clean up after run
    if [ $? -eq 0 ]
    then
        mv {rundir} {datadir}
        exit 0
    else
        mv {rundir} {faildir}
        exit 1
    fi""")


class CorsikaBatch(object):

    """Run many simultaneous CORSIKA simulations using Stoomboot

    Stoomboot is the Nikhef computer cluster.

    :param energy: the energy of the primary particle in log10(E[eV]),
                   so an energy of 16 (10**16 eV) corresponds to 10**7 GeV,
                   integer values and values ending in .5 are allowed.
    :param particle: name of primary particle, e.g. proton, gamma or iron
    :param zenith: zenith angle of the primary particle (in degrees),
                   common choices: 0, 7.5, 15, 22.5, 30, 37.5, 45 and 52.5.
    :param azimuth: azimuth angle of the primary particle (in degrees),
                    common choices: 0, 45, 90, 135, 180, 225, 270 and 315.
    :param queue: choose a queue to sumbit the job to:
                  express - max 10 minutes, max 2 jobs
                  short - max 4 hours, max 1000 jobs
                  generic - max 24 hours, max 500 jobs
                  long - max 96 hours (default 48 hours), max 500 jobs
    :param corsika: name of the compiled CORSIKA executable to use:
                    corsika74000Linux_EPOS_gheisha
                    corsika74000Linux_QGSII_gheisha
                    corsika74000Linux_QGSJET_gheisha
                    corsika74000Linux_SIBYLL_gheisha

    """

    def __init__(self, energy=16, particle='proton', zenith=22.5, azimuth=180,
                 queue='generic', corsika='corsika74000Linux_QGSII_gheisha'):
        # Energy is stored as log10(E[GeV]) for CORSIKA
        if modf(energy)[0] == 0.:
            self.energy_pre = 1.
            self.energy_pow = int(energy - 9)
        elif modf(energy)[0] == 0.5:
            self.energy_pre = 3.16228
            self.energy_pow = int(modf(energy)[1] - 9)
        else:
            raise ValueError('Energy must either be an integer or end in .5.')
        self.particle = particles.particle_id(particle)  # Store as particle id
        self.theta = zenith
        self.phi = (azimuth + 90) % 360  # Stored as Phi defined by CORSIKA
        self.queue = queue
        self.corsika = corsika
        self.seed1 = None
        self.seed2 = None
        self.rundir = None

    def run(self):
        self.prepare_env()
        self.submit_job()

    def prepare_env(self):
        """Setup CORSIKA environment"""

        # Set umask
        os.umask(002)

        # Setup directories
        taken = self.taken_seeds()
        self.generate_random_seeds(taken)
        self.make_rundir()
        self.goto_rundir()

        # Create/copy files
        self.create_input()
        self.copy_config()

    def submit_job(self):
        """Submit job to Stoomboot"""

        name = "cor_{seed1}_{seed2}".format(seed1=self.seed1, seed2=self.seed2)
        extra = "-d {rundir}".format(rundir=self.get_rundir())
        if self.queue == 'long':
            extra += " -l walltime=96:00:00"
        script = self.create_script()

        qsub.submit_job(script, name, self.queue, extra)

    def taken_seeds(self):
        """Get list of seeds already used"""

        taken = os.listdir(DATADIR)
        taken.extend(os.listdir(TEMPDIR))
        return taken

    def generate_random_seeds(self, taken):
        """Get unused combination of two seeds for CORSIKA

        :param taken: List of seed combinations already taken
                      each is formatted like this: 'seed1_seed2'

        """
        seed1 = random.randint(1, 900000000)
        seed2 = random.randint(1, 900000000)
        seed = "{seed1}_{seed2}".format(seed1=seed1, seed2=seed2)
        if seed not in taken:
            self.seed1 = seed1
            self.seed2 = seed2
            self.rundir = seed + '/'
        else:
            self.generate_random_seeds(taken)

    def make_rundir(self):
        """Make the run directory"""

        os.mkdir(self.get_rundir())

    def goto_rundir(self):
        """Move into the run directory"""

        os.chdir(self.get_rundir())

    def get_rundir(self):
        """Get run directory path"""

        return os.path.join(TEMPDIR, self.rundir)

    def create_input(self):
        """Make CORSIKA steering file"""

        input_path = os.path.join(self.get_rundir(), 'input-hisparc')
        input = INPUT_TEMPLATE.format(seed1=self.seed1, seed2=self.seed2,
                                      particle=self.particle, phi=self.phi,
                                      energy_pre=self.energy_pre,
                                      energy_pow=self.energy_pow,
                                      theta=self.theta, tablesdir=CORSIKADIR)
        with open(input_path, 'w') as input_file:
            input_file.write(input)

    def create_script(self):
        """Make Stoomboot script file"""

        exec_path = os.path.join(CORSIKADIR, self.corsika)
        run_path = self.get_rundir()

        script = SCRIPT_TEMPLATE.format(corsika=exec_path, rundir=run_path,
                                        datadir=DATADIR, faildir=FAILDIR)
        return script

    def copy_config(self):
        """Copy the CORSIKA config file to the output directory

        This way we can always check how CORSIKA was configured for each
        run.

        """
        source = os.path.join(CORSIKADIR, self.corsika + '.log')
        destination = self.get_rundir()
        subprocess.check_output(['cp', source, destination])


def multiple_jobs(n, energy, particle, zenith, azimuth, queue, corsika,
                  progress=True):
    """Use this to sumbit multiple jobs to Stoomboot

    :param n: Number of jobs to submit
    :param energy: log10(E[eV]) energy of primary particle
    :param particle: Particle kind (as string, see
                     :mod:`~sapphire.corsika.particles` for possibilities)
    :param zenith: Zenith angle in degrees of the primary particle
    :param azimuth: Azimuth angle in degrees of the primary particle
    :param queue: Stoomboot queue to submit to
    :param corsika: Name of the CORSIKA executable to use
    :param progress: Toggle printing of overview.

    """
    if progress:
        print textwrap.dedent("""\
            Batch submitting jobs to Stoomboot:
            Number of jobs      {n}
            Particle energy     10^{e} eV
            Primary particle    {p}
            Zenith angle        {z} degrees
            Azimuth angle       {a} degrees
            Stoomboot queue     {q}
            CORSIKA executable  {c}
            """.format(n=n, e=energy, p=particle, z=zenith, a=azimuth, q=queue,
                       c=corsika))

    available_slots = qsub.check_queue(queue)
    if available_slots <= 0:
        raise Exception('Submitting no jobs because selected queue is full.')
    elif available_slots < n:
        n = available_slots
        warnings.warn('Submitting {n} jobs because queue almost full.'
                      .format(n=n))

    for _ in pbar(xrange(n), show=progress):
        batch = CorsikaBatch(energy=energy, particle=particle, zenith=zenith,
                             azimuth=azimuth, queue=queue, corsika=corsika)
        batch.run()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('n', type=int, help="number of jobs to submit")
    parser.add_argument('energy', metavar='energy', type=float,
                        help="energy of the primary particle in range 11..17, "
                             "in steps of .5 (log10(E[eV]))",
                        choices=[11, 11.5, 12, 12.5, 13, 13.5, 14, 14.5,
                                 15, 15.5, 16, 16.5, 17, 17.5, 18, 18.5])
    parser.add_argument('particle', help="primary particle kind (e.g. proton "
                                         "or iron)")
    parser.add_argument('zenith', metavar='zenith',
                        help="zenith angle of primary particle in range 0..60,"
                             " in steps of 7.5 [degrees]",
                        type=float,
                        choices=[0, 7.5, 15, 22.5, 30, 37.5, 45, 52.5, 60])
    parser.add_argument('-a', '--azimuth', metavar='angle',
                        help="azimuth angle of primary particle in range "
                             "0..315, in steps of 45 [degrees]",
                        type=int,
                        default=0,
                        choices=[0, 45, 90, 135, 180, 225, 270, 315])
    parser.add_argument('-q', '--queue', metavar='name',
                        help="name of the Stoomboot queue to use, choose from "
                             "express, short, generic (default), and long",
                        default='generic',
                        choices=['express', 'short', 'generic', 'long'])
    parser.add_argument('-c', '--corsika', metavar='exec',
                        help="name of the CORSIKA executable to use",
                        default="corsika74000Linux_QGSII_gheisha")
    args = parser.parse_args()

    multiple_jobs(args.n, args.energy, args.particle, args.zenith,
                  args.azimuth, args.queue, args.corsika)


if __name__ == '__main__':
    main()
