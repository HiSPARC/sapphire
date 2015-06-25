""" Perform various angle related transformations

    Transform between different notations for angles:
    Degrees, radians and hours.

"""
from numpy import radians, degrees


def hours_to_degrees(angle):
    """Converts decimal hours to degrees

    :param hours: angle in decimal hours
    :return: angle in degrees

    """
    return angle * 15.


def hours_to_radians(angle):
    """Converts decimal hours to radians

    :param hours: angle in decimal hours
    :return: angle in radians

    """
    return radians(hours_to_degrees(angle))


def degrees_to_hours(angle):
    """Converts degrees to decimal hours

    :param angle: angle in degrees
    :return: angle in decimal hours

    """
    return angle / 15.


def radians_to_hours(angle):
    """Converts degrees to decimal hours

    :param angle: angle in degrees
    :return: angle in decimal hours

    """
    return degrees_to_hours(degrees(angle))
