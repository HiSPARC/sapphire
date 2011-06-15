import time
import hashlib
import tables

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
            print i, len(batch), hash, batch[0][0]

            data = tables.openFile(TMP_FILE % hash, 'w')
            data.createArray('/', 'positions', batch)
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

        return

    def collect_jobs(self, hashes):
        """Collect all submitted jobs and print status messages"""

        return

    def collect_results(self, hashes):
        """Collect all results into main HDF5 file"""

        return
