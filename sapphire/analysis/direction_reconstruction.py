import warnings
import itertools

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
    def reconstruct_common(cls, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0=0, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 0, 1 and 2 in ns.
        :param x# y#: position of detector 0, 1 and 2 in m.
        :param z#: height of detectors 0, 1 and 2 is ignored.

        """
        dt1 = t1 - t0
        dt2 = t2 - t0

        r1 = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        r2 = sqrt((x2 - x0) ** 2 + (y2 - x0) ** 2)

        phi1 = arctan2((y1 - y0), (x1 - x0))
        phi2 = arctan2((y2 - x0), (x2 - x0))

        return cls.reconstruct(dt1, dt2, r1, r2, phi1, phi2)

    @classmethod
    def reconstruct(cls, dt1, dt2, r1, r2, phi1, phi2):
        """Reconstruct angles from 3 detections

        :param dt#: arrival times in detector 1 and 2 relative to
                    detector 0 in ns (!).
        :param r#, phi#: position of detector 1 and 2 relative to
                         detector 0 in m and radians.
        :return: theta as given by Fokkema2012 eq 4.27,
                 phi as given by Fokkema2012 eq 4.13.

        """
        if dt1 == 0 and dt2 == 0:
            # No time difference means shower came from zenith.
            return 0, 0

        c = .3

        phi = arctan2(-(r1 * dt2 * cos(phi1) - r2 * dt1 * cos(phi2)),
                      (r1 * dt2 * sin(phi1) - r2 * dt1 * sin(phi2)))
        # The directional vector c * dt should be negative,
        # not apparent in Fokkema2012 fig 4.4.
        theta1 = arcsin(c * -dt1 / (r1 * cos(phi - phi1)))
        theta2 = arcsin(c * -dt2 / (r2 * cos(phi - phi2)))

        e1 = sqrt(cls.rel_theta1_errorsq(theta1, phi, phi1, phi2, r1, r2))
        e2 = sqrt(cls.rel_theta2_errorsq(theta2, phi, phi1, phi2, r1, r2))

        theta = (1 / e1 * theta1 + 1 / e2 * theta2) / (1 / e1 + 1 / e2)

        # We limit theta to positive values.  If theta is negative, we
        # make it positive, but need to rotate phi by 180 degrees.
        if theta < 0:
            theta *= -1
            phi += pi
            phi = (phi + pi) % (2 * pi) - pi

        if isnan(theta):
            phi = nan

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
    def reconstruct_common(cls, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0=0, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 0, 1 and 2 in ns.
        :param x# y#: position of detector 0, 1 and 2 in m.
        :param z#: height of detectors 0, 1 and 2 is ignored.

        """
        dt1 = t1 - t0
        dt2 = t2 - t0

        dx1 = x1 - x0
        dx2 = x2 - x0

        dy1 = y1 - y0
        dy2 = y2 - y0

        return cls.reconstruct(dt1, dt2, dx1, dx2, dy1, dy2)

    @staticmethod
    def reconstruct(dt1, dt2, dx1, dx2, dy1, dy2, dz1=0, dz2=0):
        """Reconstruct angles from 3 detections

        :param dt#: arrival times in detector 1 and 2 relative to
                    detector 0 in ns.
        :param dx#, dy#: position of detector 1 and 2 relative to
                         detector 0 in m.
        :param dz#: height of detectors 1 and 2 is ignored.
        :return: theta as given by Montanus2014 eq 27,
                 phi as given by Montanus2014 eq 26.

        """
        c = 0.3

        ux = c * (dt2 * dx1 - dt1 * dx2)
        uy = c * (dt2 * dy1 - dt1 * dy2)

        vz = dx1 * dy2 - dx2 * dy1

        if vz == 0:
            theta = nan
        else:
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
    def reconstruct_common(cls, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0=0, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 0, 1 and 2 in ns.
        :param x# y# z#: position of detector 0, 1 and 2 in m.

        """
        dt1 = t1 - t0
        dt2 = t2 - t0

        dx1 = x1 - x0
        dx2 = x2 - x0

        dy1 = y1 - y0
        dy2 = y2 - y0

        dz1 = z1 - z0
        dz2 = z2 - z0

        return cls.reconstruct(dt1, dt2, dx1, dx2, dy1, dy2, dz1, dz2)

    @staticmethod
    def reconstruct(dt1, dt2, dx1, dx2, dy1, dy2, dz1=0, dz2=0):
        """Reconstruct angles from 3 detections

        :param dt#: arrival times in detector 1 and 2 relative to
                    detector 0 in ns.
        :param dx#, dy#, dz#: position of detector 1 and 2 relative to
                              detector 0 in m.
        :return: theta as given by Montanus2014 eq 24,
                 phi as given by Montanus2014 eq 22.

        """
        c = .3

        ux = c * (dt2 * dx1 - dt1 * dx2)
        uy = c * (dt2 * dy1 - dt1 * dy2)
        uz = c * (dt2 * dz1 - dt1 * dz2)

        vx = dy1 * dz2 - dz1 * dy2
        vy = dz1 * dx2 - dx1 * dz2
        vz = dx1 * dy2 - dy1 * dx2

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

        if thetaplus <= pi / 2. and thetamin > pi / 2.:
            theta = thetaplus
            phi = phiplus
        elif thetaplus > pi / 2. and thetamin <= pi / 2.:
            theta = thetamin
            phi = phimin
        else:
            theta = nan
            phi = nan

        return theta, phi


class FitAlgorithm(object):

    @classmethod
    def reconstruct_common(cls, t0, t1, t2, x0, x1, x2, y0, y1, y2, z0=0, z1=0, z2=0):
        """Reconstruct angles from 3 detections

        This function converts the arguments to be suitable for the
        algorithm.

        :param t#: arrival times in detector 0, 1 and 2 in ns.
        :param x# y#: position of detector 0, 1 and 2 in m.
        :param z#: height of detectors 0, 1 and 2 is ignored.

        """
        t = (t0, t1, t2)
        x = (x0, x1, x2)
        y = (y0, y1, y2)
        z = (z0, z1, z2)

        return cls.reconstruct(t, x, y, z)

    @classmethod
    def reconstruct(cls, t, x, y, z):
        """Reconstruct angles for many detections

        :param t#: arrival times in the detectors in ns.
        :param x#, y#, z#: position of the detectors in m.
        :return: theta as given by Montanus2014 eq 21,
                 phi as given by Montanus2014 eq 22.

        """
        if not logic_checks(t, x, y, z):
            return nan, nan

        dt = cls.make_relative(t)
        dx = cls.make_relative(x)
        dy = cls.make_relative(y)
        dz = cls.make_relative(z)

        cons = ({'type': 'eq', 'fun': cls.constraint_normal_vector})

        fit = minimize(cls.best_fit, x0=(0.1, 0.1, .989),
                       args=(dt, dx, dy, dz), method="SLSQP",
                       bounds=((-1,1), (-1,1), (-1,1)), constraints=cons,
                       options={'ftol': 1e-9, 'eps': 1e-7,'maxiter': 50})
        if fit.success:
            phi1 = arctan2(fit.x[1], fit.x[0])
            theta1 = arccos(fit.x[2])
        else:
            phi1 = nan
            theta1 = nan

        fit = minimize(cls.best_fit, x0=(-0.1, -0.1, -.989),
                       args=(dt, dx, dy, dz), method="SLSQP",
                       bounds=((-1,1), (-1,1), (-1,1)), constraints=cons,
                       options={'ftol': 1e-9, 'eps': 1e-7,'maxiter': 50})
        if fit.success:
            phi2 = arctan2(fit.x[1], fit.x[0])
            theta2 = arccos(fit.x[2])
        else:
            phi2 = nan
            theta2 = nan

        # in case one of the theta's is smaller than pi/2 (shower from above)
        # and one larger than pi/2 (shower from below),
        # the first one is considered correct.
        # if both come from above (or from below), both theta's are rejected
        # the check is preceeded by a check if the fit has not delivered nans.

        if theta1 <= pi / 2. and (isnan(theta2) or theta2 > pi / 2.):
            theta = theta1
            phi = phi1
        elif (isnan(theta1) or theta1 > pi / 2.) and theta2 <= pi / 2.:
            theta = theta2
            phi = phi2
        else:
            theta = nan
            phi = nan

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
        nx, ny, nz = n_xyz

        slq = sum([(nx * xi + ny * yi + zi * nz + c * ti)**2
                   for ti, xi, yi, zi in zip(dt, dx, dy, dz)])
        return slq


class DirectEventReconstruction(DirectAlgorithmCartesian2D):

    """Reconstruct direction for station events

    This class is aware of 'events' and 'stations'.  Initialize this class
    with a 'station' and you can reconstruct events using
    :meth:`reconstruct_event`.

    :param station: :class:`sapphire.clusters.Station` object.

    """

    def __init__(self, station):
        self.station = station

    def reconstruct_event(self, event, detector_ids=[0, 2, 3]):
        """Reconstruct a single event

        :param event: an event (e.g. from an events table), or any
            dictionary-like object containing the keys necessary for
            reconstructing the direction of a shower (e.g. arrival times).
        :param detector_ids: list of the three detectors to use for
            reconstruction. The detector ids are 0-based, unlike the
            column names in the esd data.

        """
        t = [event['t%d' % (id + 1)] for id in detector_ids]
        x = [self.station.detectors[id].x for id in detector_ids]
        y = [self.station.detectors[id].y for id in detector_ids]
        z = [0, 0, 0]
        theta, phi = self.reconstruct_common(*(t + x + y + z))
        return theta, phi

    def reconstruct_events(self, events, detector_ids=[0, 2, 3]):
        """Reconstruct event

        :param event: an event (e.g. from an events table), or any
            dictionary-like object containing the keys necessary for
            reconstructing the direction of a shower (e.g. arrival times).
        :param detector_ids: detectors which use for the reconstructions.

        """
        events = events.read_where('n%d > -1 & n%d > -1 & n%d > -1' %
                                   detector_ids)
        angles = [self.reconstruct_event(event) for event in events]
        return angles


class FitClusterReconstruction(FitAlgorithm):
    def __init__(self, cluster):
        self.cluster = cluster

    def reconstruct_coincidence(coincidence):
        pass


class DirectClusterReconstruction(DirectAlgorithmCartesian3D):

    """Reconstruct coincidences with three events analytically

    This class is aware of 'coincidences' and 'clusters'.  Initialize
    this class with a 'cluster' and you can reconstruct a coincidence
    using :meth:`reconstruct_coincidence`.

    :param cluster: :class:`sapphire.clusters.Cluster` object.

    """

    def __init__(self, cluster):
        self.cluster = cluster

    def reconstruct_coincidence(coincidence, station_ids=[0, 1, 2]):
        """Reconstruct a single coincidence

        :param coincidence: a coincidence list consisting of
                            (station_number, event) tuples

        """
        pass


class ReconstructAllCoincidences():
    def __init__(self, cluster, results_table,
                 algorithm=FitClusterReconstruction):
        if algorithm == None:
            algorithm = search_algorithm()
        self.algorithm = algorithm


def logic_checks(t, x, y, z):
    """Check for impossible reconstructions

    Criteria:
    - No two detectors are at the same position.
    - Time difference between two detections should be less than distance / c.

    To add:
    - All detectors are on a line is bad.

    :param t: arrival times in the detectors in ns.
    :param x, y, z: positions of the detectors in m.
    :return: True if the checks pass, False otherwise.

    """
    # Check for identical positions
    if not len(zip(x, y, z)) == len(set(zip(x, y, z))):
        return False

    # Check if the time difference it larger than expected by c
    c = .3  # m/ns
    txyz = zip(t, x, y, z)
    for txyz0, txyz1 in itertools.combinations(txyz, 2):
        dt = abs(txyz0[0] - txyz1[0])
        dx = txyz0[1] - txyz1[1]
        dy = txyz0[2] - txyz1[2]
        dz = txyz0[3] - txyz1[3]
        dt_max = sqrt(dx**2 + dy**2 + dz**2) / c
        if dt_max < dt:
            return False

    return True
