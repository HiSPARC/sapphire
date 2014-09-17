""" Perform various axes related transformations

- Transformation between Cartesian, polar, cylindrical and spherical
  coordinate systems.
- Create a rotation matrix for rotations around a certain axis.

"""
from numpy import sqrt, arctan2, arccos, sin, cos, matrix


def cartesian_to_spherical(x, y, z):
    """Converts Cartesian coordinates into spherical coordinates

    :param x,y,z: Cartesian coordinates
    :returns: tuple of spherical coordinates (r, theta, phi),
              with theta and phi in radians.

    """
    r = sqrt(x * x + y * y + z * z)
    if r == 0:
        return 0, 0, 0
    theta = arccos(z / r)
    phi = arctan2(y, x)

    return r, theta, phi


def cartesian_to_cylindrical(x, y, z):
    """Converts Cartesian coordinates into cylindrical coordinates

    :param x,y,z: Cartesian coordinates
    :returns: tuple of cylindrical coordinates (r, phi, z), with
              phi in radians.

    """
    r = sqrt(x * x + y * y)
    phi = arctan2(y, x)
    return r, phi, z


def cartesian_to_polar(x, y):
    """Converts Cartesian coordinates into polar coordinates

    :param x,y: Cartesian coordinates
    :returns: tuple of polar coordinates (r, phi), with phi in radians.

    """
    r, phi, _ = cartesian_to_cylindrical(x, y, 0)
    return r, phi


def spherical_to_cartesian(r, theta, phi):
    """Convert spherical coordinates into Cartesian coordinates

    :param r,theta,phi: spherical coordinates, with theta and phi in radians.
    :returns: tuple of Cartesian coordinates (x, y, z)

    """
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return x, y, z


def cylindrical_to_cartesian(r, phi, z):
    """Convert cylindrical coordinates into Cartesian coordinates

    :param r,phi,z: cylindrical coordinates, with phi in radians.
    :returns: tuple of Cartesian coordinates (x, y, z)

    """
    x = r * cos(phi)
    y = r * sin(phi)
    return x, y, z


def polar_to_cartesian(r, phi):
    """Convert polar coordinates into Cartesian coordinates

    :param r,phi: polar coordinates, with phi in radians.
    :returns: tuple of Cartesian coordinates (x, y)

    """
    x, y, _ = cylindrical_to_cartesian(r, phi, 0)
    return x, y


def rotation_matrix(angle, axis='z'):
    """Generate a rotation matrix around an axis

    :param angle: amount of rotation in radians.
    :param axis: the axis to rotate around, either ``'x', 'y', 'z'``,
        or a (x,y,z) tuple specifying the axis to rotate about.

    :returns: unitary rotation matrix.

    """
    sina = sin(angle)
    cosa = cos(angle)
    if axis == 'z':
        return matrix(((cosa, sina, 0), (-sina, cosa, 0), (0, 0, 1)))
    elif axis == 'y':
        return matrix(((cosa, 0, -sina), (0, 1, 0), (sina, 0, cosa)))
    elif axis == 'x':
        return matrix(((1, 0, 0), (0, cosa, sina), (0, -sina, cosa)))
