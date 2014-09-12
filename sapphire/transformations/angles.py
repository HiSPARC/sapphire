""" Perform various angle transformations"""

from numpy import sin, cos, matrix, modf


def decimal_to_sexagesimal(decimal):
	"""Convert decimal hours or degrees to sexagesimal.

	:param decimal: decimal number to be converted to sexagismal.
	:returns: tuple of either (hours, minutes, seconds) or
	          (degrees, arcminutes, arcseconds)

	"""
	fractional, integral = modf(decimal)
	min_fractional, minutes = modf(fractional * 60)
	seconds = min_fractional * 60.
	return (integral.astype(int), minutes.astype(int), seconds)


def sexagesimal_to_decimal(hd, minutes, seconds):
	"""Convert sexagesimal hours or degrees to decimal.

	:param hd: hours or degrees.
	:param minutes: minutes or arcminutes.
	:param seconds: seconds or arcseconds.
	:returns: decimal hours or degrees.

	"""
	return hd + minutes / 60. + seconds / 3600.


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
