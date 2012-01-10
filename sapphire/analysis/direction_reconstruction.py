import sys
import os

import progressbar as pb
from numpy import isnan, arcsin, arctan2, cos, floor, inf, sin, sqrt, tan, where
from numpy.random import uniform

from sapphire import storage


class DirectionReconstruction(object):
    def __init__(self, datafile, results_table, min_n134=1., N=None, overwrite=False):
        self.data = datafile
        self.results_table = self.create_empty_output_table(results_table, overwrite)
        self.min_n134 = min_n134
        self.N = N

    def create_empty_output_table(self, table_path, overwrite=False):
        group, tablename = os.path.split(table_path)

        if not group in self.data.root:
            base, groupname = os.path.split(group)
            self.data.createGroup(base, groupname, createparents=True)
        group = self.data.getNode(group)

        if tablename in group:
            if not overwrite:
                raise RuntimeError("Reconstruction table %s already exists" % table_path)
            else:
                self.data.removeNode(group, tablename)

        table = self.data.createTable(group, tablename, storage.ReconstructedEvent)
        return table

    def reconstruct_angles_for_shower_group(self, groupname):
        """Reconstruct angles from simulation for minimum particle density"""

        shower_group = self.data.getNode(groupname)

        progressbar = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                              pb.ETA()],
                                     fd=sys.stderr)

        for shower in progressbar(self.data.listNodes(shower_group)):
            self.reconstruct_angles(shower)


    def reconstruct_angles(self, shower):
        shower_table = shower.observables
        coincidence_table = shower.coincidences
        self.station, = shower._v_attrs.cluster.stations
        if not 'cluster' in self.results_table.attrs:
            self.results_table.attrs.cluster = shower._v_attrs.cluster

        for event, coincidence in zip(shower_table[:self.N], coincidence_table[:self.N]):
            assert event['id'] == coincidence['id']
            if min(event['n1'], event['n3'], event['n4']) >= self.min_n134:
                theta, phi = self.reconstruct_angle(event)

                alpha = event['alpha']

                if not isnan(theta) and not isnan(phi):
                    self.store_reconstructed_event(coincidence, event, theta, phi)

        self.results_table.flush()

    def store_reconstructed_event(self, coincidence, event, reconstructed_theta,
                                  reconstructed_phi):
        dst_row = self.results_table.row

        dst_row['station_id'] = event['station_id']
        dst_row['r'] = coincidence['r']
        dst_row['phi'] = coincidence['phi']
        dst_row['alpha'] = event['alpha']
        dst_row['t1'] = event['t1']
        dst_row['t2'] = event['t2']
        dst_row['t3'] = event['t3']
        dst_row['t4'] = event['t4']
        dst_row['n1'] = event['n1']
        dst_row['n2'] = event['n2']
        dst_row['n3'] = event['n3']
        dst_row['n4'] = event['n4']
        dst_row['reference_theta'] = coincidence['shower_theta']
        dst_row['reference_phi'] = coincidence['shower_phi']
        dst_row['reconstructed_theta'] = reconstructed_theta
        dst_row['reconstructed_phi'] = reconstructed_phi
        dst_row['min_n134'] = min(event['n1'], event['n3'], event['n4'])
        dst_row.append()

    def reconstruct_angle(self, event):
        """Reconstruct angles from a single event"""

        c = 3.00e+8

        dt1 = event['t1'] - event['t3']
        dt2 = event['t1'] - event['t4']

        r1, phi1 = self.station.calc_r_and_phi_for_detectors(1, 3)
        r2, phi2 = self.station.calc_r_and_phi_for_detectors(1, 4)

        phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                      (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
        theta1 = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
        theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))

        e1 = sqrt(self.rel_theta1_errorsq(theta1, phi, phi1, phi2, r1, r2))
        e2 = sqrt(self.rel_theta2_errorsq(theta2, phi, phi1, phi2, r1, r2))

        theta_wgt = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

        return theta_wgt, phi

    @classmethod
    def rel_theta1_errorsq(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        sintheta = sin(theta)
        sinphiphi1 = sin(phi - phi1)

        den = (1 - sintheta ** 2) * r1 ** 2 * cos(phi - phi1) ** 2

        A = (r1 ** 2 * sinphiphi1 ** 2
             * cls.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        B = (r1 * c * sinphiphi1
             * (cls.dphi_dt0(theta, phi, phi1, phi2, r1, r2)
                - cls.dphi_dt1(theta, phi, phi1, phi2, r1, r2)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 - 2 * B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    @classmethod
    def rel_theta2_errorsq(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        sintheta = sin(theta)
        sinphiphi2 = sin(phi - phi2)

        den = (1 - sintheta ** 2) * r2 ** 2 * cos(phi - phi2) ** 2

        A = (r2 ** 2 * sinphiphi2 ** 2
             * cls.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        B = (r2 * c * sinphiphi2
             * (cls.dphi_dt0(theta, phi, phi1, phi2, r1, r2)
                - cls.dphi_dt2(theta, phi, phi1, phi2, r1, r2)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 - 2 * B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    @staticmethod
    def rel_phi_errorsq(theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) ** 2 * r1 ** 2 * r2 ** 2 * sin(theta) ** 2
           * (sinphi1 * cos(phi - phi2) - sinphi2 * cos(phi - phi1)) ** 2
           / c ** 2)

        A = (r1 ** 2 * sinphi1 ** 2
             + r2 ** 2 * sinphi2 ** 2
             - r1 * r2 * sinphi1 * sinphi2)
        B = (2 * r1 ** 2 * sinphi1 * cosphi1
             + 2 * r2 ** 2 * sinphi2 * cosphi2
             - r1 * r2 * sinphi2 * cosphi1
             - r1 * r2 * sinphi1 * cosphi2)
        C = (r1 ** 2 * cosphi1 ** 2
             + r2 ** 2 * cosphi2 ** 2
             - r1 * r2 * cosphi1 * cosphi2)

        return 2 * (A * tanphi ** 2 + B * tanphi + C) / den

    @staticmethod
    def dphi_dt0(theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta)
               * (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2))
               / c)
        num = (r2 * cosphi2 - r1 * cosphi1
               + tanphi * (r2 * sinphi2 - r1 * sinphi1))

        return num / den

    @staticmethod
    def dphi_dt1(theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta)
               * (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2))
               / c)
        num = -r2 * (sinphi2 * tanphi + cosphi2)

        return num / den

    @staticmethod
    def dphi_dt2(theta, phi, phi1, phi2, r1=10, r2=10):
        # speed of light in m / ns
        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta)
               * (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2))
               / c)
        num = r1 * (sinphi1 * tanphi + cosphi1)

        return num / den


class BinnedDirectionReconstruction(DirectionReconstruction):
    def __init__(self, datafile, results_table, min_n134=1., binning=2.5, randomize_binning=False, N=None, overwrite=False):
        super(BinnedDirectionReconstruction, self).__init__(datafile, results_table, min_n134, N, overwrite)
        self.binning = binning
        self.randomize_binning = randomize_binning

    def reconstruct_angle(self, event):
        binning = self.binning
        randomize_binning = self.randomize_binning

        for idx in 't1', 't2', 't3', 't4':
            event[idx] = floor(event[idx] / binning) * binning

        if randomize_binning is True:
            for idx in 't1', 't2', 't3', 't4':
                event[idx] += uniform(0, binning)

        return super(BinnedDirectionReconstruction, self).reconstruct_angle(event)
