import unittest

from numpy import pi, sqrt, arccos, testing, matrix

from sapphire.transformations import axes


class CoordinateSystemTests(unittest.TestCase):

    def setUp(self):
        """Test combinations of coordinates

        Cartesian, spherical, cylindrical, polar, and compass coordinates

        """
        self.combinations = (
            ((0, 0, 0), (0, 0, 0), (0, 0, 0), (0, 0, 0)),
            ((1, 0, 0), (1, pi / 2, 0), (1, 0, 0), (1, 90, 0)),
            ((-1, 0, 0), (1, pi / 2, pi), (1, pi, 0), (1, -90, 0)),
            ((0, 1, 0), (1, pi / 2, pi / 2.), (1, pi / 2., 0), (1, 0, 0)),
            ((0, -1, 0), (1, pi / 2, -pi / 2), (1, -pi / 2, 0), (1, 180, 0)),
            ((0, 0, 1), (1, 0, 0), (0, 0, 1), (0, 0, 1)),
            ((0, 0, -1), (1, pi, 0), (0, 0, -1), (0, 0, -1)),
            ((1, 1, 1), (sqrt(3), arccos(1 / sqrt(3)), pi / 4), (sqrt(2), pi / 4, 1), (sqrt(2), 45, 1)),
            ((-1, -1, -1), (sqrt(3), arccos(-1 / sqrt(3)), -pi * 3 / 4), (sqrt(2), -pi * 3 / 4, -1), (sqrt(2), -135, -1)))

    def test_cartesian_to_spherical(self):
        for cartesian, spherical, _, _ in self.combinations:
            self.assertEqual(axes.cartesian_to_spherical(*cartesian), spherical)

    def test_cartesian_to_cylindrical(self):
        for cartesian, _, cylindrical, _ in self.combinations:
            self.assertEqual(axes.cartesian_to_cylindrical(*cartesian), cylindrical)

    def test_cartesian_to_polar(self):
        for cartesian, _, cylindrical, _ in self.combinations:
            self.assertEqual(axes.cartesian_to_polar(*cartesian[:2]), cylindrical[:2])

    def test_cartesian_to_compass(self):
        for cartesian, _, _, compass in self.combinations:
            self.assertEqual(axes.cartesian_to_compass(*cartesian), compass)

    def test_spherical_to_cartesian(self):
        for cartesian, spherical, _, _ in self.combinations:
            testing.assert_almost_equal(axes.spherical_to_cartesian(*spherical), cartesian)

    def test_cylindrical_to_cartesian(self):
        for cartesian, _, cylindrical, _ in self.combinations:
            testing.assert_almost_equal(axes.cylindrical_to_cartesian(*cylindrical), cartesian)

    def test_polar_to_cartesian(self):
        for cartesian, _, cylindrical, _ in self.combinations:
            testing.assert_almost_equal(axes.polar_to_cartesian(*cylindrical[:2]), cartesian[:2])

    def test_compass_to_cartesian(self):
        for cartesian, _, _, compass in self.combinations:
            testing.assert_almost_equal(axes.compass_to_cartesian(*compass), cartesian)


class RotationMatrixTests(unittest.TestCase):

    def test_no_rotation_matrix(self):
        """Check if no rotation is correctly returned"""

        no_rotation = matrix(((1, 0, 0), (0, 1, 0), (0, 0, 1)))
        testing.assert_equal(axes.rotation_matrix(0, 'x'), no_rotation)
        testing.assert_equal(axes.rotation_matrix(0, 'y'), no_rotation)
        testing.assert_equal(axes.rotation_matrix(0, 'z'), no_rotation)

    def test_rotation_matrix(self):
        """Rotate by 90 degrees to swap the other two axes"""

        testing.assert_almost_equal(axes.rotation_matrix(pi / 2., 'x'),
                                    matrix(((1, 0, 0), (0, 0, 1), (0, -1, 0))))
        testing.assert_almost_equal(axes.rotation_matrix(pi / 2., 'y'),
                                    matrix(((0, 0, -1), (0, 1, 0), (1, 0, 0))))
        testing.assert_almost_equal(axes.rotation_matrix(pi / 2, 'z'),
                                    matrix(((0, 1, 0), (-1, 0, 0), (0, 0, 1))))


if __name__ == '__main__':
    unittest.main()
