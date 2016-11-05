"""Determine signal calibration values for data

This module can be used to determine signal calibration values from data.

Determine the PMT response curve to correct the detected number of MIPs.

"""

from numpy import log10, sqrt, where


def linear(x, slope, intercept):
    """Evaluate a linear function at point x

    :param x: x coordinate at which to evaluate the function.
    :param slope: slope of the line.
    :param intercept: y-intercept at x = 0.
    :return: the y value at x, given the slope and intercept.

    """
    return x * slope + intercept


def inverse_linear(y, slope, intercept):
    """Evaluate the inverse linear function for value y

    :param y: value y.
    :param slope: slope of the line.
    :param intercept: y-intercept at x = 0.
    :return: the x coordinate for value y, given the slope and intercept.

    """
    return (y - intercept) / slope


def log_linear(x, slope, intercept):
    """Evaluate a linear function in log-log space at point x

    :param x: x coordinate at which to evaluate the function.
    :param slope: slope of the line.
    :param intercept: log10 y-intercept at x = 1 (log10 1 = 0).
    :return: the y value at x, given the slope and intercept.

    """
    return 10 ** linear(log10(x), slope, intercept)


def linear_intersection(slope_1, intercept_1, slope_2, intercept_2):
    """The x coordinate at which the two linear lines intersect

    :param slope1: slope of the first line.
    :param intercept1: y-intercept at x = 0 of the first line.
    :param slope2: slope of the second line.
    :param intercept2: y-intercept at x = 0 of the second line.
    :return: x coordinate of the intersection.

    """
    if slope_1 == slope_2:
        raise Exception('The lines are parallel and will never cross')
    return (intercept_2 - intercept_1) / (slope_1 - slope_2)


def linear_circle_linear(x, radius, slope_low, intercept_low, slope_high, intercept_high):
    """Two linear lines connected by a circle segment

    :param x: x coordinate at which to evaluate the function.
    :param radius: radius of the cricle connecting the lines.
    :param slope_low: the slope of the lower (low x) linear line.
    :param intercept_low: the y-intercept of the lower (low x) linear line.
    :param slope_high: the slope of the upper (high x) linear line.
    :param intercept_high: the y-intercept of the upper (high x) linear line.
    :return: the y value at x.

    The center of the circle is given by (center_x, center_y)

    """
    if slope_low == slope_high:
        raise Exception('Parallel lines not allowed')

    if slope_low > slope_high:
        # Circle below the lines
        sign = -1.0
    else:
        # Circle above the lines
        sign = 1.0

    # y-intercepts of lines parallel to input lines, but shifted by r
    parallel_intercept_low = intercept_low + sign * radius * sqrt(1 + slope_low**2)
    parallel_intercept_high = intercept_high + sign * radius * sqrt(1 + slope_high**2)

    # center of the circle
    center_x = (parallel_intercept_high - parallel_intercept_low) / (slope_low - slope_high)
    center_y = slope_low * center_x + parallel_intercept_low

    # Calculate parameters for lines which intersect the circle center and
    # are perpendicual to the low or high tangents.
    perpendicular_slope_low = -1.0 / slope_low
    perpendicular_slope_high = -1.0 / slope_high
    perpendicular_intercept_low = center_x * (slope_low + 1 / slope_low) + parallel_intercept_low
    perpendicular_intercept_high = center_x * (slope_high + 1 / slope_high) + parallel_intercept_high

    # x positions where the line transition to the circle and vice-versa
    intersect_circle_low_x = linear_intersection(
        slope_low,
        intercept_low,
        perpendicular_slope_low,
        perpendicular_intercept_low,
    )
    intersect_circle_high_x = linear_intersection(
        slope_high,
        intercept_high,
        perpendicular_slope_high,
        perpendicular_intercept_high,
    )

    y = where(
        x < intersect_circle_low_x,
        linear(x, slope_low, intercept_low),
        center_y - sign * sqrt(radius**2 - (x - center_x) ** 2),
    )
    y = where(x < intersect_circle_high_x, y, linear(x, slope_high, intercept_high))

    return y


def xy_circle_linear(x, radius, slope_high, intercept_high):
    """An x=y and a linear line connected by a circle segment

    :param x: x coordinate at which to evaluate the function.
    :param radius: radius of the cricle connecting the lines.
    :param slope_high: the slope of the upper (high x) linear line.
    :param intercept_high: the y-intercept of the upper (high x) linear line.
    :return: the y value at x.

    """
    slope_low, intercept_low = (1, 0)
    return linear_circle_linear(x, radius, slope_low, intercept_low, slope_high, intercept_high)


def loglog_xy_circle_linear(x, radius, slope_high, intercept_high):
    """As `xy_circle_linear` but in log-log space"""

    return 10 ** xy_circle_linear(log10(x), radius, slope_high, intercept_high)


def loglog_linear_circle_linear(x, radius, slope_low, intercept_low, slope_high, intercept_high):
    """As `linear_circle_linear` but in log-log space"""

    return 10 ** linear_circle_linear(log10(x), radius, slope_low, intercept_low, slope_high, intercept_high)
