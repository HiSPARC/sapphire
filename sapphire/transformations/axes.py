""" Perform various axes related transformations

- Transformation between Cartesian, polar, cylindrical, spherical and compass
  coordinate systems.
- Create a rotation matrix for rotations around a certain axis.

Cartesian coordinates: x, y, z axes.

Spherical coordinates:
- r: length of vector.
- theta: angle of the vector to the z-axis.
- phi: angle of vector to the x-axis in x,y-plane, rotating counterclockwise.

Cylindrical and polar coordinates:
- r: length of vector in x,y-plane.
- phi: angle of vector to the x-axis in x,y-plane, rotating counterclockwise.
- z: height above x,y-plane.

Compass coordinates:
- r: length of vector in x,y-plane.
- alpha: angle of vector to the y-axis in x,y-plane, rotating clockwise.
- z: height above x,y-plane.

"""
from numpy import sqrt, arctan2, arccos, sin, cos, matrix, radians, degrees


def cartesian_to_spherical(x, y, z):
    """Converts Cartesian coordinates into spherical coordinates

    :param x,y,z: Cartesian coordinates
    :return: tuple of spherical coordinates (r, theta, phi),
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
    :return: tuple of cylindrical coordinates (r, phi, z), with
             phi in radians.

    """
    r = sqrt(x * x + y * y)
    phi = arctan2(y, x)
    return r, phi, z


def cartesian_to_polar(x, y):
    """Converts Cartesian coordinates into polar coordinates

    :param x,y: Cartesian coordinates
    :return: tuple of polar coordinates (r, phi), with phi in radians.

    """
    r, phi, _ = cartesian_to_cylindrical(x, y, 0)
    return r, phi


def cartesian_to_compass(x, y, z):
    """Converts Cartesian coordinates into compass coordinates

    :param x,y,z: Cartesian coordinates
    :return: tuple of compass coordinates (r, alpha, z),
             with alpha in degrees.

    """
    r = sqrt(x * x + y * y)
    alpha = degrees(arctan2(x, y))
    return r, alpha, z


def spherical_to_cartesian(r, theta, phi):
    """Convert spherical coordinates into Cartesian coordinates

    :param r,theta,phi: spherical coordinates, with theta and phi in radians.
    :return: tuple of Cartesian coordinates (x, y, z)

    """
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)
    return x, y, z


def cylindrical_to_cartesian(r, phi, z):
    """Convert cylindrical coordinates into Cartesian coordinates

    :param r,phi,z: cylindrical coordinates, with phi in radians.
    :return: tuple of Cartesian coordinates (x, y, z)

    """
    x = r * cos(phi)
    y = r * sin(phi)
    return x, y, z


def polar_to_cartesian(r, phi):
    """Convert polar coordinates into Cartesian coordinates

    :param r,phi: polar coordinates, with phi in radians.
    :return: tuple of Cartesian coordinates (x, y)

    """
    x, y, _ = cylindrical_to_cartesian(r, phi, 0)
    return x, y


def compass_to_cartesian(r, alpha, z):
    """Converts compass coordinates into Cartesian coordinates

    :param r,alpha,z: compass coordinates, with alpha in degrees.
    :return: tuple of Cartesian coordinates (x, y, z).

    """
    x = sin(radians(alpha)) * r
    y = cos(radians(alpha)) * r
    return x, y, z


def rotate_cartesian(x, y, z, angle, axis='z'):
    """Rotate Cartesian coordinates

    :param x,y,z: Cartesian coordinates.
    :param angle: amount of rotation in radians.
    :param axis: the axis to rotate around, either ``'x', 'y', 'z'``,
                 or a (x,y,z) tuple specifying the axis to rotate about.
    :return: tuple of Cartesian coordinates (x, y, z).

    """
    rot = rotation_matrix(angle, axis)
    new = (x, y, z) * rot
    return new.item(0), new.item(1), new.item(2)


def rotation_matrix(angle, axis='z'):
    """Generate a rotation matrix around an axis

    :param angle: amount of rotation in radians.
    :param axis: the axis to rotate around, either ``'x', 'y', 'z'``,
                 or a (x,y,z) tuple specifying the axis to rotate about.

    :return: unitary rotation matrix.

    """
    sina = sin(angle)
    cosa = cos(angle)
    if axis == 'z':
        return matrix(((cosa, sina, 0), (-sina, cosa, 0), (0, 0, 1)))
    elif axis == 'y':
        return matrix(((cosa, 0, -sina), (0, 1, 0), (sina, 0, cosa)))
    elif axis == 'x':
        return matrix(((1, 0, 0), (0, cosa, sina), (0, -sina, cosa)))
