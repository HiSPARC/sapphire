from numpy import (nan, isnan, arcsin, arccos, arctan2, sin, cos, tan,
                   sqrt, floor, where, deg2rad, pi, inf)


class DirectAlgorithm(object):

    """Reconstruct angles using direct analytical formula.

    This implements the equations derived in Fokkema2012 sec 4.2.
    (DOI: 10.3990/1.9789036534383)

    Note! The detectors are numbered 1-based, not 0-based.

    Speed of light is in [m / ns]

    """

    @classmethod
    def reconstruct(cls, t1, t2, t3, x1, x2, x3, y1, y2, y3, z1=0, z2=0, z3=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 1, 2 and 3 in ns.
        :param x# y#: position of detector 1, 2 and 3 in m.
        :param z#: height of detectors 1, 2 and 3 is ignored.

        """
        t2 -= t1
        t3 -= t1

        r2 = sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        r3 = sqrt((x3 - x1) ** 2 + (y3 - x1) ** 2)

        phi2 = arctan2((y2 - y1), (x2 - x1))
        phi3 = arctan2((y3 - x1), (x3 - x1))

        return cls._reconstruct(t2, t3, r2, r3, phi2, phi3)

    @classmethod
    def _reconstruct(cls, t2, t3, r2, r3, phi2, phi3):
        """Reconstruct angles from 3 detections

        :param t#: arrival times in detector 2 and 3 relative to
                   detector 1 in ns.
        :param r#, phi#: position of detector 2 and 3 relative to
                         detector 1 in m and radians.
        :return: theta as given by Fokkema2012 eq 4.27,
                 phi as given by Fokkema2012 eq 4.13.

        """
        c = .3

        t2 *= -1
        t3 *= -1

        phi = arctan2(-(r2 * t3 * cos(phi2) - r3 * t2 * cos(phi3)),
                      (r2 * t3 * sin(phi2) - r3 * t2 * sin(phi3)))
        theta2 = arcsin(c * t2 / (r2 * cos(phi - phi2)))
        theta3 = arcsin(c * t3 / (r3 * cos(phi - phi3)))

        e2 = sqrt(cls.rel_theta2_errorsq(theta2, phi, phi2, phi3, r2, r3))
        e3 = sqrt(cls.rel_theta3_errorsq(theta3, phi, phi2, phi3, r2, r3))

        theta = (1 / e2 * theta2 + 1 / e3 * theta3) / (1 / e2 + 1 / e3)

        if theta < 0:
            theta *= -1
            phi += pi
            phi = (phi + pi) % (2 * pi) - pi

        return theta, phi

    @classmethod
    def rel_theta2_errorsq(cls, theta, phi, phi2, phi3, r2=10, r3=10):
        """Fokkema2012, eq 4.23"""

        c = .3

        sintheta = sin(theta)
        sinphiphi2 = sin(phi - phi2)

        den = r2 ** 2 * (1 - sintheta ** 2) * cos(phi - phi2) ** 2

        A = (r2 ** 2 * sinphiphi2 ** 2 *
             cls.rel_phi_errorsq(theta, phi, phi2, phi3, r2, r3))
        B = -(2 * r2 * c * sinphiphi2 *
              (cls.dphi_dt1(theta, phi, phi2, phi3, r2, r3) -
               cls.dphi_dt2(theta, phi, phi2, phi3, r2, r3)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 + B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    @classmethod
    def rel_theta3_errorsq(cls, theta, phi, phi2, phi3, r2=10, r3=10):
        """Fokkema2012, eq 4.23"""

        c = .3

        sintheta = sin(theta)
        sinphiphi3 = sin(phi - phi3)

        den = r3 ** 2 * (1 - sintheta ** 2) * cos(phi - phi3) ** 2

        A = (r3 ** 2 * sinphiphi3 ** 2 *
             cls.rel_phi_errorsq(theta, phi, phi2, phi3, r2, r3))
        B = -(2 * r3 * c * sinphiphi3 *
              (cls.dphi_dt1(theta, phi, phi2, phi3, r2, r3) -
               cls.dphi_dt3(theta, phi, phi2, phi3, r2, r3)))
        C = 2 * c ** 2

        errsq = (A * sintheta ** 2 + B * sintheta + C) / den

        return where(isnan(errsq), inf, errsq)

    @staticmethod
    def rel_phi_errorsq(theta, phi, phi2, phi3, r2=10, r3=10):
        """Fokkema2012, eq 4.22"""

        c = .3

        tanphi = tan(phi)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)
        sinphi3 = sin(phi3)
        cosphi3 = cos(phi3)

        den = ((1 + tanphi ** 2) ** 2 * r2 ** 2 * r3 ** 2 * sin(theta) ** 2 *
               (sinphi2 * cos(phi - phi3) - sinphi3 * cos(phi - phi2)) ** 2 /
               c ** 2)

        A = (r2 ** 2 * sinphi2 ** 2 +
             r3 ** 2 * sinphi3 ** 2 -
             r2 * r3 * sinphi2 * sinphi3)
        B = (2 * r2 ** 2 * sinphi2 * cosphi2 +
             2 * r3 ** 2 * sinphi3 * cosphi3 -
             r2 * r3 * (sinphi3 * cosphi2 - sinphi2 * cosphi3))
        C = (r2 ** 2 * cosphi2 ** 2 +
             r3 ** 2 * cosphi3 ** 2 -
             r2 * r3 * cosphi2 * cosphi3)

        return 2 * (A * tanphi ** 2 + B * tanphi + C) / den

    @classmethod
    def dphi_dt1(cls, theta, phi, phi2, phi3, r2=10, r3=10):
        """Fokkema2012, eq 4.19"""

        return -(cls.dphi_dt2(theta, phi, phi2, phi3, r2, r3)
                 + cls.dphi_dt3(theta, phi, phi2, phi3, r2, r3))

    @staticmethod
    def dphi_dt2(theta, phi, phi2, phi3, r2=10, r3=10):
        """Fokkema2012, eq 4.20"""

        c = .3

        tanphi = tan(phi)
        sinphi2 = sin(phi2)
        sinphi3 = sin(phi3)
        cosphi3 = cos(phi3)

        den = ((1 + tanphi ** 2) * r2 * r3 * sin(theta) *
               (sinphi3 * cos(phi - phi2) - sinphi2 * cos(phi - phi3)) /
               c)
        num = -r3 * (sinphi3 * tanphi + cosphi3)

        return num / den

    @staticmethod
    def dphi_dt3(theta, phi, phi2, phi3, r2=10, r3=10):
        """Fokkema2012, eq 4.21"""

        c = .3

        tanphi = tan(phi)
        sinphi2 = sin(phi2)
        cosphi2 = cos(phi2)
        sinphi3 = sin(phi3)

        den = ((1 + tanphi ** 2) * r2 * r3 * sin(theta) *
               (sinphi3 * cos(phi - phi2) - sinphi2 * cos(phi - phi3)) /
               c)
        num = r2 * (sinphi2 * tanphi + cosphi2)

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
    def reconstruct(cls, t1, t2, t3, x1, x2, x3, y1, y2, y3, z1=0, z2=0, z3=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 1, 2 and 3 in ns.
        :param x# y#: position of detector 1, 2 and 3 in m.
        :param z#: height of detectors 1, 2 and 3 is ignored.

        """
        t2 -= t1
        t3 -= t1

        x2 -= x1
        x3 -= x1

        y2 -= y1
        y3 -= y1

        return cls._reconstruct(t2, t3, x2, x3, y2, y3)

    @staticmethod
    def _reconstruct(t2, t3, x2, x3, y2, y3, z2=0, z3=0):
        """Reconstruct angles from 3 detections

        :param t#: arrival times in detector 2 and 3 relative to
                   detector 1 in ns.
        :param x#, y#: position of detector 2 and 3 relative to
                       detector 1 in m.
        :param z#: height of detectors 2 and 3 is ignored.
        :return: theta as given by Montanus2014 eq 27,
                 phi as given by Montanus2014 eq 26.

        """
        c = 0.3

        ux = c * (t3 * x2 - t2 * x3)
        uy = c * (t3 * y2 - t2 * y3)

        vz = x2 * y3 - x3 * y2

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
    def reconstruct(cls, t1, t2, t3, x1, x2, x3, y1, y2, y3, z1=0, z2=0, z3=0):
        """Reconstruct angles from 3 detections

        This function converts the coordinates to be suitable for the
        algorithm.

        :param t#: arrival times in detector 1, 2 and 3 in ns.
        :param x# y# z#: position of detector 1, 2 and 3 in m.

        """
        t2 -= t1
        t3 -= t1

        x2 -= x1
        x3 -= x1

        y2 -= y1
        y3 -= y1

        z2 -= z1
        z3 -= z1

        return cls._reconstruct(t2, t3, x2, x3, y2, y3, z2=0, z3=0)


    @staticmethod
    def _reconstruct(t2, t3, x2, x3, y2, y3, z2=0, z3=0):
        """Reconstruct angles from 3 detections

        :param t#: arrival times in detector 2 and 3 relative to
                   detector 1 in ns.
        :param x#, y#, z#: position of detector 2 and 3 relative to
                           detector 1 in m.
        :return: theta as given by Montanus2014 eq 24,
                 phi as given by Montanus2014 eq 22.

        """
        c = .3

        ux = c * (t3 * x2 - t2 * x3)
        uy = c * (t3 * y2 - t2 * y3)
        uz = c * (t3 * z2 - t2 * z3)

        vx = y2 * z3 - y3 * z2
        vy = x3 * z2 - x2 * z3
        vz = x2 * y3 - x3 * y2

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

        phi = phiplus
        theta = thetaplus

        if thetamin < thetaplus:
            phi = phimin
            theta = thetamin

        if thetamin < pi / 2 and thetaplus < pi / 2:
            phi = nan
            theta = nan

        if thetamin > pi / 2 and thetaplus > pi / 2:
            phi = nan
            theta = nan

        return theta, phi


class FitAlgorithm(object):
    def reconstruct(times_list, positions_list):
        return theta, phi


class DirectReconstruction(DirectAlgorithm):
    def __init__(self, station):
        self.station = station

    def reconstruct_event(event):
        pass


class FitClusterReconstruction(FitAlgorithm):
    def __init__(self, cluster):
        self.cluster = cluster

    def reconstruct_coincidence(coincidence):
        pass


class ReconstructAllCoincidences():
    def __init__(self, cluster, results_table,
                 algorithm=FitClusterReconstruction):
        if algorithm == None:
            algorithm = search_algorithm()
        self.algorithm = algorithm
