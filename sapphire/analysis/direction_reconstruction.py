import warnings

from numpy import (nan, isnan, arcsin, arccos, arctan2, sin, cos, tan,
                   sqrt, floor, where, deg2rad, pi, inf)
from scipy.optimize import minimize

class DirectAlgorithm(object):

    """Reconstruct angles using direct analytical formula.

    This implements the equations derived in Fokkema2012 sec 4.2.
    (DOI: 10.3990/1.9789036534383)

    Note! The detectors are 0-based.

    Speed of light is in [m / ns]

    """

    @classmethod
    def reconstruct(cls, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0=0, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 0, 1 and 2 in ns.
        :param x# y#: position of detector 0, 1 and 2 in m.
        :param z#: height of detectors 0, 1 and 2 is ignored.

        """
        t1 -= t0
        t2 -= t0

        r1 = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        r2 = sqrt((x2 - x0) ** 2 + (y2 - x0) ** 2)

        phi1 = arctan2((y1 - y0), (x1 - x0))
        phi2 = arctan2((y2 - x0), (x2 - x0))

        return cls._reconstruct(t1, t2, r1, r2, phi1, phi2)

    @classmethod
    def _reconstruct(cls, t1, t2, r1, r2, phi1, phi2):
        """Reconstruct angles from 3 detections

        :param t#: arrival times in detector 1 and 2 relative to
                   detector 0 in ns.
        :param r#, phi#: position of detector 1 and 2 relative to
                         detector 0 in m and radians.
        :return: theta as given by Fokkema2012 eq 4.27,
                 phi as given by Fokkema2012 eq 4.13.

        """
        c = .3

        t1 = -t1
        t2 = -t2

        phi = arctan2(-(r1 * t2 * cos(phi1) - r2 * t1 * cos(phi2)),
                      (r1 * t2 * sin(phi1) - r2 * t1 * sin(phi2)))
        theta1 = arcsin(c * t1 / (r1 * cos(phi - phi1)))
        theta2 = arcsin(c * t2 / (r2 * cos(phi - phi2)))

        e1 = sqrt(cls.rel_theta1_errorsq(theta1, phi, phi1, phi2, r1, r2))
        e2 = sqrt(cls.rel_theta2_errorsq(theta2, phi, phi1, phi2, r1, r2))

        theta = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

        if theta < 0:
            theta *= -1
            phi += pi
            phi = (phi + pi) % (2 * pi) - pi

        return theta, phi

    @classmethod
    def rel_theta1_errorsq(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.23"""

        c = .3

        sintheta = sin(theta)
        sinphiphi1 = sin(phi - phi1)

        den = r1 ** 2 * (1 - sintheta ** 2) * cos(phi - phi1) ** 2

        A = (r1 ** 2 * sinphiphi1 ** 2 *
             cls.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        B = -(2 * r1 * c * sinphiphi1 *
              (cls.dphi_dt0(theta, phi, phi1, phi2, r1, r2) -
               cls.dphi_dt1(theta, phi, phi1, phi2, r1, r2)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 + B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    @classmethod
    def rel_theta2_errorsq(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.23"""

        c = .3

        sintheta = sin(theta)
        sinphiphi2 = sin(phi - phi2)

        den = r2 ** 2 * (1 - sintheta ** 2) * cos(phi - phi2) ** 2

        A = (r2 ** 2 * sinphiphi2 ** 2 *
             cls.rel_phi_errorsq(theta, phi, phi1, phi2, r1, r2))
        B = -(2 * r2 * c * sinphiphi2 *
              (cls.dphi_dt0(theta, phi, phi1, phi2, r1, r2) -
               cls.dphi_dt2(theta, phi, phi1, phi2, r1, r2)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 + B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    @staticmethod
    def rel_phi_errorsq(theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.22"""

        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) ** 2 * r1 ** 2 * r2 ** 2 * sin(theta) ** 2 *
               (sinphi1 * cos(phi - phi2) - sinphi2 * cos(phi - phi1)) ** 2 /
               c ** 2)

        A = (r1 ** 2 * sinphi1 ** 2 +
             r2 ** 2 * sinphi2 ** 2 -
             r1 * r2 * sinphi1 * sinphi2)
        B = (2 * r1 ** 2 * sinphi1 * cosphi1 +
             2 * r2 ** 2 * sinphi2 * cosphi2 -
             r1 * r2 * (sinphi2 * cosphi1 - sinphi1 * cosphi2))
        C = (r1 ** 2 * cosphi1 ** 2 +
             r2 ** 2 * cosphi2 ** 2 -
             r1 * r2 * cosphi1 * cosphi2)

        return 2 * (A * tanphi ** 2 + B * tanphi + C) / den

    @classmethod
    def dphi_dt0(cls, theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.19"""

        return -(cls.dphi_dt1(theta, phi, phi1, phi2, r1, r2) +
                 cls.dphi_dt2(theta, phi, phi1, phi2, r1, r2))

    @staticmethod
    def dphi_dt1(theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.20"""

        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta) *
               (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2)) /
               c)
        num = -r2 * (sinphi2 * tanphi + cosphi2)

        return num / den

    @staticmethod
    def dphi_dt2(theta, phi, phi1, phi2, r1=10, r2=10):
        """Fokkema2012, eq 4.21"""

        c = .3

        tanphi = tan(phi)
        sinphi1 = sin(phi1)
        cosphi1 = cos(phi1)
        sinphi2 = sin(phi2)

        den = ((1 + tanphi ** 2) * r1 * r2 * sin(theta) *
               (sinphi2 * cos(phi - phi1) - sinphi1 * cos(phi - phi2)) /
               c)
        num = r1 * (sinphi1 * tanphi + cosphi1)

        return num / den


class DirectAlgorithmCartesian2D(object):

    """Reconstruct angles using direct analytical formula.

    This implements the equations derived in Montanus2014.
    "Direction reconstruction of cosmic air showers with
     detectorstations at different altitudes"

    Here the 2D version is used, assuming each detector is at the same
    altitude.

    """

    @classmethod
    def reconstruct(cls, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0=0, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 0, 1 and 2 in ns.
        :param x# y#: position of detector 0, 1 and 2 in m.
        :param z#: height of detectors 0, 1 and 2 is ignored.

        """
        t1 -= t0
        t2 -= t0

        x1 -= x0
        x2 -= x0

        y1 -= y0
        y2 -= y0

        return cls._reconstruct(t1, t2, x1, x2, y1, y2)

    @staticmethod
    def _reconstruct(t1, t2, x1, x2, y1, y2, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        :param t#: arrival times in detector 1 and 2 relative to
                   detector 0 in ns.
        :param x#, y#: position of detector 1 and 2 relative to
                       detector 0 in m.
        :param z#: height of detectors 1 and 2 is ignored.
        :return: theta as given by Montanus2014 eq 27,
                 phi as given by Montanus2014 eq 26.

        """
        c = 0.3

        ux = c * (t2 * x1 - t1 * x2)
        uy = c * (t2 * y1 - t1 * y2)

        vz = x1 * y2 - x2 * y1

        usquared = ux * ux + uy * uy
        vzsquared = vz * vz

        phi = arctan2(-ux * vz, uy * vz)
        theta = arcsin(sqrt(usquared / vzsquared))

        if isnan(theta):
            phi = nan

        return theta, phi


class DirectAlgorithmCartesian3D(object):

    """Reconstruct angles using direct analytical formula.

    This implements the equations derived in Montanus2014.
    "Direction reconstruction of cosmic air showers with
     detectorstations at different altitudes"

    Here the 3D version is used, assuming each detector is at the same
    altitude.

    """

    @classmethod
    def reconstruct(cls, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0=0, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 0, 1 and 2 in ns.
        :param x# y# z#: position of detector 0, 1 and 2 in m.

        """
        t1 -= t0
        t2 -= t0

        x1 -= x0
        x2 -= x0

        y1 -= y0
        y2 -= y0

        z1 -= z0
        z2 -= z0

        return cls._reconstruct(t1, t2, x1, x2, y1, y2, z1=0, z2=0)


    @staticmethod
    def _reconstruct(t1, t2, x1, x2, y1, y2, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        :param t#: arrival times in detector 1 and 2 relative to
                   detector 0 in ns.
        :param x#, y#, z#: position of detector 1 and 2 relative to
                           detector 0 in m.
        :return: theta as given by Montanus2014 eq 24,
                 phi as given by Montanus2014 eq 22.

        """
        c = .3

        ux = c * (t2 * x1 - t1 * x2)
        uy = c * (t2 * y1 - t1 * y2)
        uz = c * (t2 * z1 - t1 * z2)

        vx = y1 * z2 - z1 * y2
        vy = z1 * x2 - x1 * z2
        vz = x1 * y2 - y1 * x2

        ucrossvx = uy * vz - uz * vy
        ucrossvy = uz * vx - ux * vz
        ucrossvz = ux * vy - uy * vx

        usquared = ux * ux + uy * uy + uz * uz
        vsquared = vx * vx + vy * vy + vz * vz
        underroot = vsquared - usquared
        if underroot < 0:
            underroot = nan

        termx = vx * sqrt(underroot)
        termy = vy * sqrt(underroot)
        termz = vz * sqrt(underroot)

        nxplus = (ucrossvx + termx) / vsquared
        nyplus = (ucrossvy + termy) / vsquared
        nzplus = (ucrossvz + termz) / vsquared

        nxmin = (ucrossvx - termx) / vsquared
        nymin = (ucrossvy - termy) / vsquared
        nzmin = (ucrossvz - termz) / vsquared

        phiplus = arctan2(nyplus, nxplus)
        thetaplus = arccos(nzplus)

        phimin = arctan2(nymin, nxmin)
        thetamin = arccos(nzmin)

        if isnan(thetaplus):
            thetaplus = pi

        if isnan(thetamin):
            thetamin = pi

        if thetamin < thetaplus:
            phi = phimin
            theta = thetamin
        else:
            phi = phiplus
            theta = thetaplus

        if thetamin < pi / 2 and thetaplus < pi / 2:
            phi = nan
            theta = nan
        elif thetamin > pi / 2 and thetaplus > pi / 2:
            phi = nan
            theta = nan

        return theta, phi


class FitAlgorithm(object):

    @classmethod
    def reconstruct(cls, t, x, y, z):
        """Reconstruct angles for many detections

        :param t#: arrival times in the detectors in ns.
        :param x#, y#, z#: position of the detectors in m.
        :return: theta as given by Montanus2014 eq 21,
                 phi as given by Montanus2014 eq 22.

        """
        dt = cls.make_relative(t)
        dx = cls.make_relative(x)
        dy = cls.make_relative(y)
        dz = cls.make_relative(z)

        cons = ({'type': 'eq', 'fun': cls.constraint_normal_vector})
        fit = minimize(cls.best_fit, x0=(0, 0, 1), args=(dt, dx, dy, dz),
                       method="SLSQP", bounds=((-1,1), (-1,1), (-1,1)),
                       constraints=cons)
        phi = arctan2(fit.x[1], fit.x[0])
        theta = arccos(fit.x[2])

        if theta < 0:
            warnings.warn('Theta was less than 0')
            theta *= -1
            phi += pi
            phi = (phi + pi) % (2 * pi) - pi

        return theta, phi

    @staticmethod
    def make_relative(x):
        """Make first element the origin and make rest relative to it."""

        return [xi - x[0] for xi in x[1:]]

    @staticmethod
    def constraint_normal_vector(n):
        """This should be equal to zero"""

        return n[0]**2 + n[1]**2 + n[2]**2 - 1

    @staticmethod
    def best_fit(n_xyz, dt, dx, dy, dz):
        """The function to be minimized to find the direction

        :param n_xyz: list containing the unit vector.
        :param dt: list of relative arrival times in the detectors in ns.
        :param dx, dy, dz: list of relative detector positions in m.
        :return: least sum of squares as in Montanus2014, eq 36

        """
        c = .3
        nx = n_xyz[0]
        ny = n_xyz[1]
        nz = n_xyz[2]

        slq = sum([(nx * xi + ny * yi + zi * nz + c * ti)**2
                   for ti, xi, yi, zi in zip(dt, dx, dy, dz)])
        return slq


class DirectReconstruction(DirectAlgorithm):

    """Reconstruct event using :class:`DirectAlgorithm`

    This class is aware of 'events' and 'stations'.  Initialize this class
    with a 'station' and you can reconstruct events using
    :meth:`reconstruct_event`.

    :param station: :class:`sapphire.clusters.Station` object.

    """

    def __init__(self, station):
        self.station = station

    def reconstruct_event(event, detector_ids=[0, 2, 3]):
        """Reconstruct a single event

        :param event: an event (e.g. from an events table), or any
            dictionary-like object containing the keys necessary for
            reconstructing the direction of a shower (e.g. arrival times).

        """
        pass


class FitClusterReconstruction(FitAlgorithm):
    def __init__(self, cluster):
        self.cluster = cluster

    def reconstruct_coincidence(coincidence):
        pass


class DirectClusterReconstruction(DirectAlgorithm):

    """Reconstruct coincidence using :class:`DirectAlgorithm`

    This class is aware of 'events' and 'clusters'.  Initialize this class
    with a 'cluster' and you can reconstruct coincidences using
    :meth:`reconstruct_coincidence`.

    :param cluster: :class:`sapphire.clusters.Cluster` object.

    """

    def __init__(self, cluster):
        self.cluster = cluster

    def reconstruct_coincidence(coincidence, station_ids=[0, 1, 2]):
        """Reconstruct a single coincidence

        :param coincidence: a coincidence (e.g. from a coincidences
            table)

        """
        pass


class ReconstructAllCoincidences():
    def __init__(self, cluster, results_table,
                 algorithm=FitClusterReconstruction):
        if algorithm == None:
            algorithm = search_algorithm()
        self.algorithm = algorithm
