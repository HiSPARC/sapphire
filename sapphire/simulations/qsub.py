import time
import hashlib
import tables
import os
import subprocess
import progressbar as pb

from groundparticles import GroundParticlesSimulation
from sapphire import storage

JOB_FILE = '__QSUB_%s.h5'
STATUS_FILE = '__STATUS_%s'
STAT_INTERVAL = 2
BAR_LENGTH = 40


class QSubSimulation(GroundParticlesSimulation):

    """Submit a simulation in multiple jobs

    This class splits up the simulation in several jobs.  Each job is
    submitted through qsub and the results are collected.
    """

    def __init__(self, cluster, data, grdpcles, output, R, N, N_cores=2,
                 *args, **kwargs):
        """Initialize the simulation

        :param cluster: BaseCluster (or derived) instance
        :param data: the HDF5 file
        :param grdpcles: name of the dataset containing the ground particles
        :param output: name of the destination group to store results
        :param R: maximum distance of shower to center of cluster
        :param N: number of simulations to perform
        :param N_cores: number of cores to run on

        """
        self.N_cores = N_cores
        super(QSubSimulation, self).__init__(cluster, data, grdpcles,
                                             output, R, N, *args, **kwargs)

    def run(self):
        """Perform a simulation

        Perform a simulation by initializing and submitting a set of qsub
        jobs.  When all jobs are done, collect all data and store it.

        Note that the simulation is split into N_cores jobs if, and only
        if, N divided by N_cores is an integer.  If not, the simulation is
        split into N_cores + 1 jobs, with the last job containing the
        remainder of events.

        """
        positions = list(self.generate_positions())

        hashes = []
        for i, batch in enumerate(self._chunks(positions)):
            hash = hashlib.sha1(str(time.time()) + str(i)).hexdigest()[:8]
            hashes.append(hash)

            data = tables.open_file(JOB_FILE % hash, 'w')
            data.create_array('/', 'positions', batch)
            data.root._v_attrs.cluster = self.cluster
            data.root._v_attrs.data = self.data.filename
            data.root._v_attrs.grdpcles = self.grdpcles._v_pathname
            data.root._v_attrs.output = self.output._v_pathname
            data.root._v_attrs.N = len(batch)
            data.root._v_attrs.R = self.R
            data.root._v_attrs.use_poisson = self.use_poisson
            data.root._v_attrs.gauss = self.gauss
            data.root._v_attrs.trig_threshold = self.trig_threshold
            data.close()
            self._qsub(hash)

        self.collect_jobs(hashes)
        self.collect_results(hashes)

        self.store_observables()

    def _chunks(self, data):
        """Split data into chunks max N_cores + 1 chunks"""

        N = int(self.N / self.N_cores)
        for i in xrange(0, len(data), N):
            yield data[i:i + N]

    def _qsub(self, hash):
        """Submit a batch job using qsub"""

        dir = os.path.dirname(__file__)
        environ = os.environ.copy()
        environ['JOB_HASH'] = hash
        try:
            output = subprocess.check_output(['qsub', '-V', '-q', 'short',
                                              '%s/qsub.sh' % dir],
                                              env=environ)
        except subprocess.CalledProcessError, exc:
            print 80 * '-'
            print exc
            print "Program output given below:"
            print exc.output
        else:
            print output

    def collect_jobs(self, hashes):
        """Collect all submitted jobs and print status messages"""

        N_T = len(hashes)
        progressbar = pb.ProgressBar(maxval=100 * N_T,
                                     widgets=[pb.Percentage(), pb.Bar(),
                                              pb.ETA()]).start()
        old_progress = 0

        while True:
            N_Q = 0
            N_R = 0
            N_C = 0
            msgs = []
            progress = 0

            for hash in hashes:
                status = STATUS_FILE % hash
                if os.path.exists(status):
                    with open(status, 'r') as f:
                        status = f.read()
                    if status.find('DONE') >= 0:
                        # job finished
                        N_C += 1
                        progress += 100
                    else:
                        # job running
                        N_R += 1

                        try:
                            status = int(status)
                            progress += status
                        except ValueError:
                            continue

                        bar = int(BAR_LENGTH * status / 100)
                        bar = bar * '#' + (BAR_LENGTH - bar) * ' '
                        msgs.append('job %s: %2d %% |%s|' % (hash, status,
                                                             bar))
                else:
                    # job queued
                    N_Q += 1

            # ANSI: ED (Erase Data, or: clear screen)
            print chr(27) + "[2J"
            # ANSI: CUP (CUrsor Position)
            print chr(27) + "[;H"

            print "Status Display"
            print
            if progress > old_progress:
                progressbar.update(progress)
                old_progress = progress
            print
            print
            print "Queued: %d/%d, Running: %d/%d, Completed: %d/%d" % \
                (N_Q, N_T, N_R, N_T, N_C, N_T)
            print 80 * '-'
            for msg in msgs:
                print msg

            if N_C == N_T:
                # all jobs finished
                print "Simulation finished!"
                progressbar.finish()
                print
                return

            time.sleep(2)

    def collect_results(self, hashes):
        """Collect all results into main HDF5 file"""

        headers = self.data.create_table(self.output, '_headers',
                                         storage.SimulationEventHeader)
        particles = self.data.create_table(self.output, '_particles',
                                         storage.SimulationParticle)

        base_id = 0
        for hash in hashes:
            job_data = tables.open_file(JOB_FILE % hash, 'r')
            job_headers = job_data.get_node(self.output, '_headers').read()
            job_particles = job_data.get_node(self.output,
                                             '_particles').read()

            job_headers['id'] += base_id
            job_particles['id'] += base_id

            headers.append(job_headers)
            particles.append(job_particles)

            base_id += job_data.root._v_attrs.N
            job_data.close()
            os.remove(JOB_FILE % hash)
            os.remove(STATUS_FILE % hash)

        self.headers = headers
        self.particles = particles


class QSubChild(GroundParticlesSimulation):

    """Accept a job from QSubSimulation and run it

    This class accepts a job submitted through QSubSimulation and runs it.
    """

    def __init__(self, data, hash):
        """Initialize the simulation

        This is a bit different from a regular simulation in that we have
        one shower data file (read-only) and several job-specific data
        files.  The job-specific data file (identified by the hash)
        contains several simulation parameters as attributes on the root
        group.

        :param data: job HDF5 file
        :param hash: job hash

        """
        self.data = data
        self.hash = hash

        attrs = data.root._v_attrs
        self.cluster = attrs.cluster
        self.N = attrs.N
        self.R = attrs.R
        self.use_poisson = attrs.use_poisson
        self.gauss = attrs.gauss
        self.trig_threshold = attrs.trig_threshold

        output = attrs.output
        head, tail = os.path.split(output)
        self.output = self.data.create_group(head, tail,
                                            createparents=True)

        self.attr_shower_data = attrs.data
        self.attr_grdpcles = attrs.grdpcles

    def run(self):
        """Run the simulation

        Fetch the positions to be simulated from the data file and run the
        simulation

        """
        shower_data = tables.open_file(self.attr_shower_data, 'r')
        self.grdpcles = shower_data.get_node(self.attr_grdpcles)

        self._run_welcome_msg()
        positions = self.data.root.positions.read()

        self.headers = self.data.create_table(self.output, '_headers',
                                             storage.SimulationEventHeader)
        self.particles = self.data.create_table(self.output, '_particles',
                                               storage.SimulationParticle)

        N = 0
        prev = time.time()
        for event_id, (r, phi, alpha) in enumerate(positions):
            self.simulate_event(event_id, r, phi, alpha)
            N += 1
            now = time.time()
            if now - prev > STAT_INTERVAL:
                prev = now
                with open(STATUS_FILE % hash, 'w') as f:
                    f.write("%d" % (100 * N / self.N))

        self.headers.flush()
        self.particles.flush()
        shower_data.close()

        self._run_exit_msg()


if __name__ == '__main__':
    hash = os.environ['JOB_HASH']
    print "Running job", hash

    data = tables.open_file(JOB_FILE % hash, 'a')
    sim = QSubChild(data, hash)
    sim.run()
    data.close()

    with open(STATUS_FILE % hash, 'w') as f:
        f.write('DONE')
