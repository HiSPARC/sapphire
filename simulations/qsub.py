import time
import hashlib
import tables
import os
import subprocess

from base import BaseSimulation

TMP_FILE = '__QSUB_%s.h5'


class QSubSimulation(BaseSimulation):

    """Submit a simulation in multiple jobs

    This class splits up the simulation in several jobs.  Each job is
    submitted through qsub and the results are collected.
    """

    def __init__(self, cluster, data, grdpcles, output, R, N, N_cores=2):
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
                                             output, R, N)

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

            data = tables.openFile(TMP_FILE % hash, 'w')
            data.createArray('/', 'positions', batch)
            data.root._v_attrs.cluster = self.cluster
            data.root._v_attrs.data = self.data.filename
            data.root._v_attrs.grdpcles = self.grdpcles._v_pathname
            data.root._v_attrs.output = self.output._v_pathname
            data.root._v_attrs.N = len(batch)
            data.root._v_attrs.R = self.R
            data.close()
            self._qsub(hash)

        self.collect_jobs(hashes)
        self.collect_results(hashes)

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
            output = subprocess.check_output(['%s/qsub.sh' % dir],
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

        return

    def collect_results(self, hashes):
        """Collect all results into main HDF5 file"""

        return


class QSubChild(BaseSimulation):

    """Accept a job from QSubSimulation and run it

    This class accepts a job submitted through QSubSimulation and runs it.
    """

    def __init__(self, hash):
        """Initialize the simulation

        This is a bit different from a regular simulation in that we have
        one shower data file (read-only) and several job-specific data
        files.  The job-specific data file (identified by the hash)
        contains several simulation parameters as attributes on the root
        group.

        :param hash: job hash

        """
        data = tables.openFile(TMP_FILE % hash, 'a')
        attrs = data.root._v_attrs

        self.data = data
        self.cluster = attrs.cluster
        self.N = attrs.N
        self.R = attrs.R

        output = attrs.output
        head, tail = os.path.split(output)
        self.output = self.data.createGroup(head, tail,
                                            createparents=True)

        shower_data = tables.openFile(attrs.data, 'r')
        self.shower_data = shower_data
        self.grdpcles = shower_data.getNode(attrs.grdpcles)

    def run(self):
        """Run the simulation

        Fetch the positions to be simulated from the data file and run the
        simulation

        """
        positions = self.data.root.positions.read()
        self._do_run(positions)


if __name__ == '__main__':
    hash = os.environ['JOB_HASH']
    print "Running job", hash

    sim = QSubChild(hash)
    sim.run()
