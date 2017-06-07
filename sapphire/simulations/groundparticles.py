"""Perform simulations of CORSIKA air showers on a cluster of stations

This simulation uses a HDF5 file created from a CORSIKA simulation with
the ``store_corsika_data`` script. The shower is 'thrown' on the cluster
with random core positions and azimuth angles.

Example usage::

    >>> import tables
    >>> from sapphire import GroundParticlesSimulation, ScienceParkCluster
    >>> data = tables.open_file('/tmp/test_groundparticle_simulation.h5', 'w')
    >>> cluster = ScienceParkCluster()
    >>> sim = GroundParticlesSimulation('corsika.h5', 500, cluster, data,
    ...                                 '/', 10)
    >>> sim.run()

"""
from __future__ import print_function

from math import pi, sin, cos, tan, sqrt, log10
from time import time
import subprocess
import shutil
import os
from six import iteritems

import numpy as np
import tables

from .gammas import simulate_detector_mips_gammas
from .detector import HiSPARCSimulation, ErrorlessSimulation
from ..corsika.corsika_queries import CorsikaQuery
from ..corsika.particles import particle_id
from ..utils import pbar, norm_angle, closest_in_list, vector_length, c


class GroundParticlesGEANT4Simulation(HiSPARCSimulation):

    def __init__(self, corsikafile_path, max_core_distance, *args, **kwargs):
        """Simulation initialization

        :param corsikafile_path: path to the corsika.h5 file containing
                                 the groundparticles.
        :param max_core_distance: maximum distance of shower core to
                                  center of cluster.

        """
        super(GroundParticlesGEANT4Simulation, self).__init__(*args, **kwargs)

        self.corsikafile = tables.open_file(corsikafile_path, 'r')
        self.groundparticles = self.corsikafile.get_node('/groundparticles')
        self.max_core_distance = max_core_distance

    def __del__(self):
        self.finish()

    def finish(self):
        """Clean-up after simulation"""

        self.corsikafile.close()

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple) and azimuth.

        """
        r_max = self.max_core_distance
        now = int(time())

        event_header = self.corsikafile.get_node_attr('/', 'event_header')
        event_end = self.corsikafile.get_node_attr('/', 'event_end')
        corsika_parameters = {'zenith': event_header.zenith,
                              'size': event_end.n_electrons_levels,
                              'energy': event_header.energy,
                              'particle': event_header.particle}
        self.corsika_azimuth = event_header.azimuth

        self.corsika_zenith = corsika_parameters['zenith']
        self.corsika_energy = corsika_parameters['energy']
        self.cr_particle = particle_id(corsika_parameters['particle'])

        for i in pbar(range(self.n), show=self.progress):
            ext_timestamp = (now + i) * int(1e9)
            x, y = self.generate_core_position(r_max)
            self.core_distance = np.sqrt(x**2 + y**2)
            self.shower_azimuth = self.generate_azimuth()

            shower_parameters = {'ext_timestamp': ext_timestamp,
                                 'core_pos': (x, y),
                                 'azimuth': self.shower_azimuth}

            # Subtract CORSIKA shower azimuth from desired shower azimuth
            # make it fit in (-pi, pi] to get rotation angle of the cluster.
            alpha = self.shower_azimuth - self.corsika_azimuth
            alpha = norm_angle(alpha)
            self._prepare_cluster_for_shower(x, y, alpha)

            shower_parameters.update(corsika_parameters)
            yield shower_parameters

    def _prepare_cluster_for_shower(self, x, y, alpha):
        """Prepare the cluster object for the simulation of a shower.

        Rotate and translate the cluster so that (0, 0) coincides with the
        shower core position and that the angle between the rotated cluster
        and the CORSIKA shower is the desired azimuth.

        :param x,y: position of shower core relative to cluster origin in m.
        :param alpha: angle the cluster needs to be rotated in radians.

        """
        # rotate the core position around the original cluster center
        xp = x * cos(-alpha) - y * sin(-alpha)
        yp = x * sin(-alpha) + y * cos(-alpha)

        self.cluster.set_coordinates(-xp, -yp, 0, -alpha)

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Checks if leptons have passed a detector. If so, it returns the number
        of leptons in the detector and the arrival time of the first lepton
        passing the detector.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         the observables will be determined.
        :param shower_parameters: dictionary with the shower parameters.

        """

        particles = self.get_particles_in_detector(detector, shower_parameters)
        n_detected = len(particles)

        if n_detected:
            n_muons, n_electrons, n_gammas, firstarrival, pulseintegral, \
            pulseintegral_muon, pulseintegral_electron, pulseintegral_gamma = \
            self.simulate_detector_mips_for_particles(particles, detector, 
                                                      shower_parameters)
            particles['t'] += firstarrival
            nz = cos(shower_parameters['zenith'])
            tproj = detector.get_coordinates()[-1] / (c * nz)
            first_signal = particles['t'].min() + detector.offset - tproj
            observables = {'n': n_muons + n_electrons + n_gammas,
                           'n_muons': n_muons,
                           'n_electrons': n_electrons,
                           'n_gammas': n_gammas,
                           't': self.simulate_adc_sampling(first_signal),
                           'integrals': pulseintegral,
                           'integrals_muon': pulseintegral_muon,
                           'integrals_electron': pulseintegral_electron,
                           'integrals_gamma': pulseintegral_gamma}
        else:
            observables = {'n': 0, 'n_muons': 0, 'n_electrons': 0, 'n_gammas': 0,
                           't': -999, 'integrals': 0., 'integrals_muon': 0.,
                           'integrals_electron': 0., 'integrals_gamma': 0.}

        return observables

    def simulate_detector_mips_for_particles(self, particles, detector,
                                             shower_parameters):
        """Simulate the detector signal for particles

        :param particles: particle rows with the p_[x, y, z]
                          components of the particle momenta.

        """
        
        # Run the geant4 simulation for each particle
        arrived_photons_per_particle = []
        arrived_photons_per_particle_muon = []
        arrived_photons_per_particle_electron = []
        arrived_photons_per_particle_gamma = []
        arrivaltimes = []
        n_muons = 0
        n_electrons = 0
        n_gammas = 0
        particle_types = ['', 'gamma', 'e+', 'e-', '', 'mu+', 'mu-']
        for particle in particles:
            # Determine which particle hit the detector
            particle_id = particle["particle_id"]
            particletype = particle_types[particle_id]
            
            # Determine the position the particle hit the detector in the
            # detector reference system (-25 < x < 25 and -50 < y < 50)
            # taking projection due to detector-height differences into
            # account.
            x = particle["x"]
            y = particle["y"]
            p = np.array([x, y])
            
            detx, dety, detz = detector.get_coordinates()
            detcorners = detector.get_corners()
            
            zenith = shower_parameters['zenith']
            azimuth = self.corsika_azimuth

            znxnz = detz * tan(zenith) * cos(azimuth)
            znynz = detz * tan(zenith) * sin(azimuth)
            
            detcproj = [(cx - znxnz, cy - znynz) for cx, cy in detcorners]
            cproj1 = np.array([detcproj[0][0],detcproj[0][1]])
            cproj2 = np.array([detcproj[1][0],detcproj[1][1]])
            cproj3 = np.array([detcproj[2][0],detcproj[2][1]])
            cproj4 = np.array([detcproj[3][0],detcproj[3][1]])
            
            # Here I determine the distance from a point to a line
            ydistance = (np.linalg.norm(np.cross(cproj2 - cproj1, cproj1 - p)) /
                         np.linalg.norm(cproj2 - cproj1))
            xdistance = (np.linalg.norm(np.cross(cproj4 - cproj1, cproj1 - p)) /
                         np.linalg.norm(cproj4 - cproj1))
            
            xdetcoord = 100 * xdistance - 25
            ydetcoord = 100 * ydistance - 50
            
            # Determine at which angle the particle hit the detector
            px = particle["p_x"]
            py = particle["p_y"]
            pz = particle["p_z"]

            # Determine the energy of the incoming particle
            particleenergy = np.sqrt(px ** 2 + py ** 2 + pz ** 2)
            
            # Start the GEANT4 simulation using the position, direction and
            # energy of the incoming particle. This simulation creates a
            # new directory RUN_1 with a csv file containing the number of
            # photons that arrived at the PMT.
            output = subprocess.check_output(["/user/kaspervd/Documents/repositories/diamond/20170117_geant4_simulation/HiSPARC-stbc-build/./skibox", "1", particletype,
                                              "{}".format(particleenergy),
                                              "{}".format(ydetcoord),
                                              "{}".format(xdetcoord),
                                              "-99899",#"-99893.695",
                                              "{}".format(px),
                                              "{}".format(py),
                                              "{}".format(pz)])
            #"""
            print( "./skibox", "1", particletype,
                                              "{}".format(particleenergy),
                                              "{}".format(ydetcoord),
                                              "{}".format(xdetcoord),
                                              "-99893.695",
                                              "{}".format(px),
                                              "{}".format(py),
                                              "{}".format(pz) )
            #"""

            # Determine the number of photons that have arrived at the PMT
            # and the time it took for the first photon to arrive at the PMT.
            file = np.genfromtxt("RUN_1/outpSD.csv", delimiter=",")
            try:
                numberofphotons = len(file[1:, 1]) # first element is header
                arrivaltime = min(file[1:, 0])
                # Succesful interaction, keep statistics
                if particle_id == 1:
                    print("Gamma detected, strength: {}".format(numberofphotons))
                    n_gammas += 1
                    arrived_photons_per_particle_gamma.append(numberofphotons)
                elif particle_id in [2, 3]:
                    print("Electron detected, strength: {}".format(numberofphotons))
                    n_electrons += 1
                    arrived_photons_per_particle_electron.append(numberofphotons)
                elif particle_id in [5, 6]:
                    print("Muon detected, strength: {}".format(numberofphotons))
                    n_muons += 1
                    arrived_photons_per_particle_muon.append(numberofphotons)
            except:
                # No photons have arrived (a gamma that didn't undergo any
                # iteraction).
                numberofphotons = 0
                arrivaltime = -999
            
            # Remove the directory created by the GEANT4 simulation
            shutil.rmtree("RUN_1")

            # If multiple particles hit the detector, they are treated
            # seperately. Make lists in order to be able to add all
            # arrived photons.
            arrived_photons_per_particle.append(numberofphotons)
            arrivaltimes.append(arrivaltime)

        # Combine the separate geant4 runs to obtain the signal
        # (in arrived photons)
        pulseintegral = np.sum(np.array(arrived_photons_per_particle))
        pulseintegral_muon = np.sum(np.array(arrived_photons_per_particle_muon))
        pulseintegral_electron = np.sum(np.array(arrived_photons_per_particle_electron))
        pulseintegral_gamma = np.sum(np.array(arrived_photons_per_particle_gamma))
        firstarrival = min(arrivaltimes)
    
        # Mip here stands for arrived photons.
        return n_muons, n_electrons, n_gammas, firstarrival, pulseintegral, \
               pulseintegral_muon, pulseintegral_electron, pulseintegral_gamma
    
    def simulate_trigger(self, detector_observables):
        """Simulate a trigger response.

        This implements the trigger as used on HiSPARC stations:
        - 4-detector station: at least two high or three low signals.
        - 2-detector station: at least 2 low signals.

        :param detector_observables: list of dictionaries, each containing
                                     the observables of one detector.
        :return: True if the station triggers, False otherwise.

        """
        n_detectors = len(detector_observables)
        detectors_low = sum([True for observables in detector_observables
                             if observables['integrals'] > 30])
        detectors_high = sum([True for observables in detector_observables
                              if observables['integrals'] > 70])

        if n_detectors == 4 and (detectors_high >= 2 or detectors_low >= 3):
            return True
        elif n_detectors == 2 and detectors_low >= 2:
            return True
        else:
            return False

    def simulate_gps(self, station_observables, shower_parameters, station):
        """Simulate gps timestamp.

        :param station_observables: dictionary containing the observables
                                    of the station.
        :param shower_parameters: dictionary with the shower parameters.
        :param station: :class:`sapphire.clusters.Station` for which
                         to simulate the gps timestamp.
        :return: station_observables updated with gps timestamp and
                 trigger time.

        """
        arrival_times = [station_observables['t%d' % id]
                         for id in range(1, 5)
                         if station_observables.get('n%d' % id, -1) > 0]

        if len(arrival_times) > 1:
            trigger_time = sorted(arrival_times)[1]

            ext_timestamp = shower_parameters['ext_timestamp']
            ext_timestamp += int(trigger_time + station.gps_offset +
                                 self.simulate_gps_uncertainty())
            timestamp = int(ext_timestamp / int(1e9))
            nanoseconds = int(ext_timestamp % int(1e9))

            gps_timestamp = {'ext_timestamp': ext_timestamp,
                             'timestamp': timestamp,
                             'nanoseconds': nanoseconds,
                             't_trigger': trigger_time}
            station_observables.update(gps_timestamp)

        return station_observables

    def get_particles_in_detector(self, detector, shower_parameters):
        """Simulate the detector detection area accurately.

        First particles are filtered to see which fall inside a
        non-rotated square box around the detector (i.e. sides of 1.2m).
        For the remaining particles a more accurate query is used to see
        which actually hit the detector. The advantage of using the
        square is that column indexes can be used, which may speed up
        queries.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.
        :param shower_parameters: dictionary with the shower parameters.

        """
        
        # Possible keys in particles
        #
        # particle_id, 1 = gamma, 2-3 is electron, 4 is neutrino, 5-6 is muon
        # r - core distance in m
        # phi - azimuth angle in rad
        # x - x position in m
        # y - y position in m
        # t - time since first interaction in ns
        # p_x - momentum in x direction in eV/c
        # p_y - momentum in y direction in eV/c
        # p_z - momentum in z direction in eV/c
        # hadron_generation
        # observation_level - observation level above sea level in cm
        
        detector_boundary = 0.6

        x, y, z = detector.get_coordinates()
        corners = detector.get_corners()
        zenith = shower_parameters['zenith']
        azimuth = self.corsika_azimuth

        znxnz = z * tan(zenith) * cos(azimuth)
        znynz = z * tan(zenith) * sin(azimuth)
        xproj = x - znxnz
        yproj = y - znynz

        cproj = [(cx - znxnz, cy - znynz) for cx, cy in corners]

        b11, line1, b12 = self.get_line_boundary_eqs(*cproj[0:3])
        b21, line2, b22 = self.get_line_boundary_eqs(*cproj[1:4])
        query = ("(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f) & "
                 "(b11 < %s) & (%s < b12) & (b21 < %s) & (%s < b22) & "
                 "(particle_id <= 6)" %
                 (xproj - detector_boundary, xproj + detector_boundary,
                  yproj - detector_boundary, yproj + detector_boundary,
                  line1, line1, line2, line2))

        return self.groundparticles.read_where(query)

    def get_line_boundary_eqs(self, p0, p1, p2):
        """Get line equations using three points

        Given three points, this function computes the equations for two
        parallel lines going through these points.  The first and second
        point are on the same line, whereas the third point is taken to
        be on a line which runs parallel to the first.  The return value
        is an equation and two boundaries which can be used to test if a
        point is between the two lines.

        :param p0,p1: (x, y) tuples on the same line.
        :param p2: (x, y) tuple on the parallel line.
        :return: value1, equation, value2, such that points satisfying
            value1 < equation < value2 are between the parallel lines.

        Example::

            >>> get_line_boundary_eqs((0, 0), (1, 1), (0, 2))
            (0.0, 'y - 1.000000 * x', 2.0)

        """
        (x0, y0), (x1, y1), (x2, y2) = p0, p1, p2

        # Compute the general equation for the lines
        if x0 == x1:
            # line is exactly vertical
            line = "x"
            b1, b2 = x0, x2
        else:
            # First, compute the slope
            a = (y1 - y0) / (x1 - x0)

            # Calculate the y-intercepts of both lines
            b1 = y0 - a * x0
            b2 = y2 - a * x2

            line = "y - %f * x" % a

        # And order the y-intercepts
        if b1 > b2:
            b1, b2 = b2, b1

        return b1, line, b2

    def store_station_observables(self, station_id, station_observables):
        """Store station observables.

        :param station_id: the id of the station in self.cluster
        :param station_observables: A dictionary containing the
            variables to be stored for this event.
        :return: The index (row number) of the newly added event.

        """
        events_table = self.station_groups[station_id].events
        row = events_table.row
        row['event_id'] = events_table.nrows
        row['shower_energy'] = self.corsika_energy
        row['zenith'] = self.corsika_zenith
        row['azimuth'] = self.shower_azimuth
        row['cr_particle'] = self.cr_particle
        row['core_distance'] = self.core_distance
        for key, value in iteritems(station_observables):
            if key in events_table.colnames:
                row[key] = value
            else:
                warnings.warn('Unsupported variable')
        row.append()
        events_table.flush()

        return events_table.nrows - 1


class GroundParticlesSimulation(HiSPARCSimulation):

    def __init__(self, corsikafile_path, max_core_distance, *args, **kwargs):
        """Simulation initialization

        :param corsikafile_path: path to the corsika.h5 file containing
                                 the groundparticles.
        :param max_core_distance: maximum distance of shower core to
                                  center of cluster.

        """
        super(GroundParticlesSimulation, self).__init__(*args, **kwargs)

        self.corsikafile = tables.open_file(corsikafile_path, 'r')
        self.groundparticles = self.corsikafile.get_node('/groundparticles')
        self.max_core_distance = max_core_distance

    def __del__(self):
        self.finish()

    def finish(self):
        """Clean-up after simulation"""

        self.corsikafile.close()

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple) and azimuth.

        """
        r_max = self.max_core_distance
        now = int(time())

        event_header = self.corsikafile.get_node_attr('/', 'event_header')
        event_end = self.corsikafile.get_node_attr('/', 'event_end')
        corsika_parameters = {'zenith': event_header.zenith,
                              'size': event_end.n_electrons_levels,
                              'energy': event_header.energy,
                              'particle': event_header.particle}
        self.corsika_azimuth = event_header.azimuth

        for i in pbar(range(self.n), show=self.progress):
            ext_timestamp = (now + i) * int(1e9)
            x, y = self.generate_core_position(r_max)
            shower_azimuth = self.generate_azimuth()

            shower_parameters = {'ext_timestamp': ext_timestamp,
                                 'core_pos': (x, y),
                                 'azimuth': shower_azimuth}

            # Subtract CORSIKA shower azimuth from desired shower azimuth
            # make it fit in (-pi, pi] to get rotation angle of the cluster.
            alpha = shower_azimuth - self.corsika_azimuth
            alpha = norm_angle(alpha)
            self._prepare_cluster_for_shower(x, y, alpha)

            shower_parameters.update(corsika_parameters)
            yield shower_parameters

    def _prepare_cluster_for_shower(self, x, y, alpha):
        """Prepare the cluster object for the simulation of a shower.

        Rotate and translate the cluster so that (0, 0) coincides with the
        shower core position and that the angle between the rotated cluster
        and the CORSIKA shower is the desired azimuth.

        :param x,y: position of shower core relative to cluster origin in m.
        :param alpha: angle the cluster needs to be rotated in radians.

        """
        # rotate the core position around the original cluster center
        xp = x * cos(-alpha) - y * sin(-alpha)
        yp = x * sin(-alpha) + y * cos(-alpha)

        self.cluster.set_coordinates(-xp, -yp, 0, -alpha)

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Checks if leptons have passed a detector. If so, it returns the number
        of leptons in the detector and the arrival time of the first lepton
        passing the detector.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         the observables will be determined.
        :param shower_parameters: dictionary with the shower parameters.

        """

        particles = self.get_particles_in_detector(detector, shower_parameters)
        n_detected = len(particles)

        if n_detected:
            mips = self.simulate_detector_mips_for_particles(particles)
            particles['t'] += self.simulate_signal_transport_time(n_detected)
            nz = cos(shower_parameters['zenith'])
            tproj = detector.get_coordinates()[-1] / (c * nz)
            first_signal = particles['t'].min() + detector.offset - tproj
            observables = {'n': round(mips, 3),
                           't': self.simulate_adc_sampling(first_signal)}
        else:
            observables = {'n': 0., 't': -999}

        return observables

    def simulate_detector_mips_for_particles(self, particles):
        """Simulate the detector signal for particles

        :param particles: particle rows with the p_[x, y, z]
                          components of the particle momenta.

        """
        # determination of lepton angle of incidence
        theta = np.arccos(abs(particles['p_z']) /
                          vector_length(particles['p_x'], particles['p_y'],
                                        particles['p_z']))
        n = len(particles)
        mips = self.simulate_detector_mips(n, theta)

        return mips

    def simulate_trigger(self, detector_observables):
        """Simulate a trigger response.

        This implements the trigger as used on HiSPARC stations:
        - 4-detector station: at least two high or three low signals.
        - 2-detector station: at least 2 low signals.

        :param detector_observables: list of dictionaries, each containing
                                     the observables of one detector.
        :return: True if the station triggers, False otherwise.

        """
        n_detectors = len(detector_observables)
        detectors_low = sum([True for observables in detector_observables
                             if observables['n'] > 0.3])
        detectors_high = sum([True for observables in detector_observables
                              if observables['n'] > 0.5])

        if n_detectors == 4 and (detectors_high >= 2 or detectors_low >= 3):
            return True
        elif n_detectors == 2 and detectors_low >= 2:
            return True
        else:
            return False

    def simulate_gps(self, station_observables, shower_parameters, station):
        """Simulate gps timestamp.

        :param station_observables: dictionary containing the observables
                                    of the station.
        :param shower_parameters: dictionary with the shower parameters.
        :param station: :class:`sapphire.clusters.Station` for which
                         to simulate the gps timestamp.
        :return: station_observables updated with gps timestamp and
                 trigger time.

        """
        arrival_times = [station_observables['t%d' % id]
                         for id in range(1, 5)
                         if station_observables.get('n%d' % id, -1) > 0]

        if len(arrival_times) > 1:
            trigger_time = sorted(arrival_times)[1]

            ext_timestamp = shower_parameters['ext_timestamp']
            ext_timestamp += int(trigger_time + station.gps_offset +
                                 self.simulate_gps_uncertainty())
            timestamp = int(ext_timestamp / int(1e9))
            nanoseconds = int(ext_timestamp % int(1e9))

            gps_timestamp = {'ext_timestamp': ext_timestamp,
                             'timestamp': timestamp,
                             'nanoseconds': nanoseconds,
                             't_trigger': trigger_time}
            station_observables.update(gps_timestamp)

        return station_observables

    def get_particles_in_detector(self, detector, shower_parameters):
        """Get particles that hit a detector.

        Particle ids 2, 3, 5, 6 are electrons and muons,
        id 4 is no longer used (were neutrino's).

        The detector is approximated by a square with a surface of 0.5
        square meter which is *not* correctly rotated.  In fact, during
        the simulation, the rotation of the detector is undefined.  This
        is faster than a more thorough implementation.

        The CORSIKA simulation azimuth is used for the projection because the
        cluster is rotated such that from the perspective of the rotated
        detectors the CORSIKA showers come from the desired azimuth. In the
        simulation frame the CORSIKA shower azimuth remains unchanged.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.
        :param shower_parameters: dictionary with the shower parameters.

        """
        detector_boundary = sqrt(0.5) / 2.

        x, y, z = detector.get_coordinates()
        zenith = shower_parameters['zenith']
        azimuth = self.corsika_azimuth

        nxnz = tan(zenith) * cos(azimuth)
        nynz = tan(zenith) * sin(azimuth)
        xproj = x - z * nxnz
        yproj = y - z * nynz

        query = ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
                 ' & (particle_id >= 2) & (particle_id <= 6)' %
                 (xproj - detector_boundary, xproj + detector_boundary,
                  yproj - detector_boundary, yproj + detector_boundary))
        return self.groundparticles.read_where(query)


class GroundParticlesGammaSimulation(GroundParticlesSimulation):
    """Simulation which includes signals from gamma particles in the shower"""

    def simulate_detector_response(self, detector, shower_parameters):
        """Simulate detector response to a shower.

        Checks if particles have passed a detector. If so, it returns the
        number of particles in the detector and the arrival time of the first
        particle passing the detector.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         the observables will be determined.
        :param shower_parameters: dictionary with the shower parameters.

        """
        leptons, gammas = self.get_particles_in_detector(detector,
                                                         shower_parameters)
        n_leptons = len(leptons)
        n_gammas = len(gammas)

        if not n_leptons + n_gammas:
            return {'n': 0, 't': -999}

        if n_leptons:
            mips_lepton = self.simulate_detector_mips_for_particles(leptons)
            leptons['t'] += self.simulate_signal_transport_time(n_leptons)
            first_lepton = leptons['t'].min()
        else:
            mips_lepton = 0

        if n_gammas:
            mips_gamma = self.simulate_detector_mips_for_gammas(gammas)
            gammas['t'] += self.simulate_signal_transport_time(n_gammas)
            first_gamma = gammas['t'].min()
        else:
            mips_gamma = 0

        if n_leptons and n_gammas:
            first_signal = min(first_lepton, first_gamma) + detector.offset
        elif n_leptons:
            first_signal = first_lepton + detector.offset
        elif n_gammas:
            first_signal = first_gamma + detector.offset

        return {'n': mips_lepton + mips_gamma,
                't': self.simulate_adc_sampling(first_signal)}

    def get_particles_in_detector(self, detector, shower_parameters):
        """Get particles that hit a detector.

        Particle ids 2, 3, 5, 6 are electrons and muons,
        id 4 is no longer used (were neutrino's).

        The detector is approximated by a square with a surface of 0.5
        square meter which is *not* correctly rotated.  In fact, during
        the simulation, the rotation of the detector is undefined.  This
        is faster than a more thorough implementation.

        *Detector height is ignored!*

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.
        :param shower_parameters: dictionary with the shower parameters.

        """
        detector_boundary = sqrt(.5) / 2.

        x, y, z = detector.get_coordinates()
        zenith = shower_parameters['zenith']
        azimuth = self.corsika_azimuth

        nxnz = tan(zenith) * cos(azimuth)
        nynz = tan(zenith) * sin(azimuth)
        xproj = x - z * nxnz
        yproj = y - z * nynz

        query_leptons = \
            ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
             ' & (particle_id >= 2) & (particle_id <= 6)' %
             (xproj - detector_boundary, xproj + detector_boundary,
              yproj - detector_boundary, yproj + detector_boundary))

        query_gammas = \
            ('(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f)'
             ' & (particle_id == 1)' %
             (xproj - detector_boundary, xproj + detector_boundary,
              yproj - detector_boundary, yproj + detector_boundary))

        return (self.groundparticles.read_where(query_leptons),
                self.groundparticles.read_where(query_gammas))

    def simulate_detector_mips_for_gammas(self, particles):
        """Simulate the detector signal for gammas

        :param particles: particle rows with the p_[x, y, z]
                          components of the particle momenta.

        """
        p_gamma = np.sqrt(particles['p_x'] ** 2 + particles['p_y'] ** 2 +
                          particles['p_z'] ** 2)

        # determination of lepton angle of incidence
        theta = np.arccos(abs(particles['p_z']) /
                          p_gamma)

        mips = simulate_detector_mips_gammas(p_gamma, theta)

        return mips


class DetectorBoundarySimulation(GroundParticlesSimulation):

    """More accuratly simulate the detection area of the detectors.

    Take the orientation of the detectors into account and use the
    exact detector boundaries. This requires a slightly more complex
    query which is a bit slower.

    """

    def get_particles_in_detector(self, detector, shower_parameters):
        """Simulate the detector detection area accurately.

        First particles are filtered to see which fall inside a
        non-rotated square box around the detector (i.e. sides of 1.2m).
        For the remaining particles a more accurate query is used to see
        which actually hit the detector. The advantage of using the
        square is that column indexes can be used, which may speed up
        queries.

        :param detector: :class:`~sapphire.clusters.Detector` for which
                         to get particles.
        :param shower_parameters: dictionary with the shower parameters.

        """
        detector_boundary = 0.6

        x, y, z = detector.get_coordinates()
        corners = detector.get_corners()
        zenith = shower_parameters['zenith']
        azimuth = self.corsika_azimuth

        znxnz = z * tan(zenith) * cos(azimuth)
        znynz = z * tan(zenith) * sin(azimuth)
        xproj = x - znxnz
        yproj = y - znynz

        cproj = [(cx - znxnz, cy - znynz) for cx, cy in corners]

        b11, line1, b12 = self.get_line_boundary_eqs(*cproj[0:3])
        b21, line2, b22 = self.get_line_boundary_eqs(*cproj[1:4])
        query = ("(x >= %f) & (x <= %f) & (y >= %f) & (y <= %f) & "
                 "(b11 < %s) & (%s < b12) & (b21 < %s) & (%s < b22) & "
                 "(particle_id >= 2) & (particle_id <= 6)" %
                 (xproj - detector_boundary, xproj + detector_boundary,
                  yproj - detector_boundary, yproj + detector_boundary,
                  line1, line1, line2, line2))

        return self.groundparticles.read_where(query)

    def get_line_boundary_eqs(self, p0, p1, p2):
        """Get line equations using three points

        Given three points, this function computes the equations for two
        parallel lines going through these points.  The first and second
        point are on the same line, whereas the third point is taken to
        be on a line which runs parallel to the first.  The return value
        is an equation and two boundaries which can be used to test if a
        point is between the two lines.

        :param p0,p1: (x, y) tuples on the same line.
        :param p2: (x, y) tuple on the parallel line.
        :return: value1, equation, value2, such that points satisfying
            value1 < equation < value2 are between the parallel lines.

        Example::

            >>> get_line_boundary_eqs((0, 0), (1, 1), (0, 2))
            (0.0, 'y - 1.000000 * x', 2.0)

        """
        (x0, y0), (x1, y1), (x2, y2) = p0, p1, p2

        # Compute the general equation for the lines
        if x0 == x1:
            # line is exactly vertical
            line = "x"
            b1, b2 = x0, x2
        else:
            # First, compute the slope
            a = (y1 - y0) / (x1 - x0)

            # Calculate the y-intercepts of both lines
            b1 = y0 - a * x0
            b2 = y2 - a * x2

            line = "y - %f * x" % a

        # And order the y-intercepts
        if b1 > b2:
            b1, b2 = b2, b1

        return b1, line, b2


class ParticleCounterSimulation(GroundParticlesSimulation):

    """Do not simulate mips, just count the number of particles."""

    def simulate_detector_mips(self, n, theta):
        """A mip for a mip, count number of particles in a detector."""

        return n


class FixedCoreDistanceSimulation(GroundParticlesSimulation):

    """Shower core at a fixed core distance (from cluster origin).

    :param core_distance: distance of shower core to center of cluster.

    """

    @classmethod
    def generate_core_position(cls, r_max):
        """Generate a random core position on a circle

        :param r_max: Fixed core distance, in meters.
        :return: Random x, y position on the circle with radius r_max.

        """
        phi = np.random.uniform(-pi, pi)
        x = r_max * cos(phi)
        y = r_max * sin(phi)
        return x, y


class GroundParticlesSimulationWithoutErrors(ErrorlessSimulation,
                                             GroundParticlesSimulation):

    """This simulation does not simulate errors/uncertainties

    This results in perfect timing (first particle through detector)
    and particle counting for the detectors.

    """

    pass


class MultipleGroundParticlesSimulation(GroundParticlesSimulation):

    """Use multiple CORSIKA simulated air showers in one run.

    Simulations will be selected from the set of available showers.
    Each time an energy and zenith angle is generated a shower is selected
    from the CORSIKA overview. Each shower is reused multiple times to
    take advantage of caching, and to reduce IO stress.

    .. warning::

        This simulation loads a new shower often it is therefore more I/O
        intensive than :class:`GroundParticlesSimulation`. Do not run many
        of these simulations simultaneously!

    """

    # CORSIKA data location at Nikhef
    DATA = '/data/hisparc/corsika/data/{seeds}/corsika.h5'

    def __init__(self, corsikaoverview_path, max_core_distance, min_energy,
                 max_energy, *args, **kwargs):
        """Simulation initialization

        :param corsikaoverview_path: path to the corsika_overview.h5 file
                                     containing the available simulations.
        :param max_core_distance: maximum distance of shower core to
                                  center of cluster.
        :param min_energy,max_energy: upper and lower shower energy limits,
                                      in eV.

        """
        # Super of the super class.
        super(GroundParticlesSimulation, self).__init__(*args, **kwargs)

        self.cq = CorsikaQuery(corsikaoverview_path)
        self.max_core_distance = max_core_distance
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.available_energies = {e for e in self.cq.all_energies
                                   if min_energy <= 10 ** e <= max_energy}
        self.available_zeniths = {e: self.cq.available_parameters('zenith',
                                                                  energy=e)
                                  for e in self.available_energies}

    def finish(self):
        """Clean-up after simulation"""

        self.cq.finish()

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple) and azimuth.

        """
        r = self.max_core_distance
        n_reuse = 100
        now = int(time())

        for i in pbar(range(self.n), show=self.progress):
            sim = self.select_simulation()
            if sim is None:
                continue

            corsika_parameters = {'zenith': sim['zenith'],
                                  'size': sim['n_electron'],
                                  'energy': sim['energy'],
                                  'particle': sim['particle_id']}
            self.corsika_azimuth = sim['azimuth']

            seeds = self.cq.seeds([sim])[0]
            with tables.open_file(self.DATA.format(seeds=seeds), 'r') as data:
                try:
                    self.groundparticles = data.get_node('/groundparticles')
                except tables.NoSuchNodeError:
                    print('No groundparticles in %s' % seeds)
                    continue

                for j in range(n_reuse):
                    ext_timestamp = (now + i + (float(j) / n_reuse)) * int(1e9)
                    x, y = self.generate_core_position(r)
                    shower_azimuth = self.generate_azimuth()

                    shower_parameters = {'ext_timestamp': ext_timestamp,
                                         'core_pos': (x, y),
                                         'azimuth': shower_azimuth}

                    # Subtract CORSIKA shower azimuth from desired shower
                    # azimuth to get rotation angle of the cluster.
                    alpha = shower_azimuth - self.corsika_azimuth
                    alpha = norm_angle(alpha)
                    self._prepare_cluster_for_shower(x, y, alpha)

                    shower_parameters.update(corsika_parameters)
                    yield shower_parameters

    def select_simulation(self):
        """Generate parameters for selecting a CORSIKA simulation

        :return: simulation row from a CORSIKA Simulations table.

        """
        energy = self.generate_energy(self.min_energy, self.max_energy)
        shower_energy = closest_in_list(log10(energy), self.available_energies)

        zenith = self.generate_zenith()
        shower_zenith = closest_in_list(np.degrees(zenith),
                                        self.available_zeniths[shower_energy])

        sims = self.cq.simulations(energy=shower_energy, zenith=shower_zenith)
        if not len(sims):
            return None
        sim = np.random.choice(sims)
        return sim


class MultipleGroundParticlesGEANT4Simulation(GroundParticlesGEANT4Simulation):

    """Use multiple CORSIKA simulated air showers in one run.

    Simulations will be selected from the set of available showers.
    Each time an energy and zenith angle is generated a shower is selected
    from the CORSIKA overview. Each shower is reused multiple times to
    take advantage of caching, and to reduce IO stress.

    .. warning::

        This simulation loads a new shower often it is therefore more I/O
        intensive than :class:`GroundParticlesSimulation`. Do not run many
        of these simulations simultaneously!

    """

    # CORSIKA data location at Nikhef
    DATA = '/dcache/hisparc/corsika/data/{seeds}/corsika.h5'

    def __init__(self, corsikaoverview_path, max_core_distance, min_energy,
                 max_energy, *args, **kwargs):
        """Simulation initialization

        :param corsikaoverview_path: path to the corsika_overview.h5 file
                                     containing the available simulations.
        :param max_core_distance: maximum distance of shower core to
                                  center of cluster.
        :param min_energy,max_energy: upper and lower shower energy limits,
                                      in eV.

        """
        super(GroundParticlesGEANT4Simulation, self).__init__(*args, **kwargs)

        self.cq = CorsikaQuery(corsikaoverview_path)
        self.max_core_distance = max_core_distance
        self.min_energy = min_energy
        self.max_energy = max_energy
        self.available_energies = {e for e in self.cq.all_energies
                                   if min_energy <= 10 ** e <= max_energy}
        self.available_zeniths = {e: self.cq.available_parameters('zenith',
                                                                  energy=e)
                                  for e in self.available_energies}

    def finish(self):
        """Clean-up after simulation"""

        self.cq.finish()

    def generate_shower_parameters(self):
        """Generate shower parameters like core position, energy, etc.

        For this groundparticles simulation, only the shower core position
        and rotation angle of the shower are generated.  Do *not*
        interpret these parameters as the position of the cluster, or the
        rotation of the cluster!  Interpret them as *shower* parameters.

        :return: dictionary with shower parameters: core_pos
                 (x, y-tuple) and azimuth.

        """
        r = self.max_core_distance
        n_reuse = 1
        now = int(time())

        for i in range(self.n):
            sim = self.select_simulation()
            if sim is None:
                continue

            corsika_parameters = {'zenith': sim['zenith'],
                                  'size': sim['n_electron'],
                                  'energy': sim['energy'],
                                  'particle': sim['particle_id']}
            self.corsika_azimuth = sim['azimuth']

            self.corsika_zenith = sim['zenith']
            self.corsika_energy = sim['energy']
            self.cr_particle = sim['particle_id']

            seeds = self.cq.seeds([sim])[0]

            if self.corsika_energy < (9.9*10**14):
                # Because of the high dcache i/o load I create, all
                # CORSIKA simulations with an energy below log(eV) = 15 
                # were moved to a temporary directory on the stoomboot node.
                tmpdir = os.environ["TMPDIR"]
                localDATA = tmpdir+"/data/{seeds}/corsika.h5"
                #print("Load local")

                with tables.open_file(localDATA.format(seeds=seeds), 'r') as data:
                    try:
                        self.groundparticles = data.get_node('/groundparticles')
                    except tables.NoSuchNodeError:
                        print('No groundparticles in %s' % seeds)
                        continue

                    for j in range(n_reuse):
                        ext_timestamp = (now + i + (float(j) / n_reuse)) * int(1e9)
                        x, y = self.generate_core_position(r)
                        self.core_distance = np.sqrt(x**2 + y**2)
                        self.shower_azimuth = self.generate_azimuth()

                        shower_parameters = {'ext_timestamp': ext_timestamp,
                                             'core_pos': (x, y),
                                             'azimuth': self.shower_azimuth}

                        # Subtract CORSIKA shower azimuth from desired shower
                        # azimuth to get rotation angle of the cluster.
                        alpha = self.shower_azimuth - self.corsika_azimuth
                        alpha = norm_angle(alpha)
                        self._prepare_cluster_for_shower(x, y, alpha)
    
                        shower_parameters.update(corsika_parameters)
                        yield shower_parameters

            else: # Use the regular dcache data
                #print("Load dCache")
                with tables.open_file(self.DATA.format(seeds=seeds), 'r') as data:
                    try:
                        self.groundparticles = data.get_node('/groundparticles')
                    except tables.NoSuchNodeError:
                        print('No groundparticles in %s' % seeds)
                        continue

                    for j in range(n_reuse):
                        ext_timestamp = (now + i + (float(j) / n_reuse)) * int(1e9)
                        x, y = self.generate_core_position(r)
                        self.core_distance = np.sqrt(x**2 + y**2)
                        self.shower_azimuth = self.generate_azimuth()

                        shower_parameters = {'ext_timestamp': ext_timestamp,
                                             'core_pos': (x, y),
                                             'azimuth': self.shower_azimuth}

                        # Subtract CORSIKA shower azimuth from desired shower
                        # azimuth to get rotation angle of the cluster.
                        alpha = self.shower_azimuth - self.corsika_azimuth
                        alpha = norm_angle(alpha)
                        self._prepare_cluster_for_shower(x, y, alpha)
    
                        shower_parameters.update(corsika_parameters)
                        yield shower_parameters


    def select_simulation(self):
        """Generate parameters for selecting a CORSIKA simulation

        :return: simulation row from a CORSIKA Simulations table.

        """
        energy = self.generate_energy(self.min_energy, self.max_energy)
        shower_energy = closest_in_list(log10(energy), self.available_energies)

        zenith = self.generate_zenith()
        shower_zenith = closest_in_list(np.degrees(zenith),
                                        self.available_zeniths[shower_energy])

        sims = self.cq.simulations(energy=shower_energy, zenith=shower_zenith)
        if not len(sims):
            return None
        sim = np.random.choice(sims)
        return sim
