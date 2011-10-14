from math import pi, sqrt
import unittest

from mock import Mock, patch

import clusters


class DetectorTests(unittest.TestCase):
    def setUp(self):
        self.mock_station = Mock()
        self.detector_1 = clusters.Detector(self.mock_station, 1, 0, 'LR')
        self.detector_2 = clusters.Detector(self.mock_station, -1, 2, 'UD')

    def test_detector_size(self):
        self.assertEqual(self.detector_1.detector_size, (.5, 1.))
        self.assertEqual(self.detector_2.detector_size, (.5, 1.))

    def test_attributes(self):
        self.assertIs(self.detector_1.station, self.mock_station)
        self.assertEqual(self.detector_1.x, 1)
        self.assertEqual(self.detector_1.y, 0)
        self.assertEqual(self.detector_1.orientation, 'LR')

    def test_LR_get_corners(self):
        self.mock_station.get_xyalpha_coordinates.return_value = (.25, 3, 0)
        corners = self.detector_1.get_corners()
        self.assertEqual(corners, [(.75, 2.75), (1.75, 2.75), (1.75, 3.25), (.75, 3.25)])

    def test_LR_get_corners_rotated(self):
        self.mock_station.get_xyalpha_coordinates.return_value = (0, 0, pi / 2)
        corners = self.detector_1.get_corners()
        expected_corners = [(.25, .5), (.25, 1.5), (-.25, 1.5), (-.25, .5)]
        for (x, y), (expected_x, expected_y) in zip(corners, expected_corners):
            self.assertAlmostEqual(x, expected_x)
            self.assertAlmostEqual(y, expected_y)

    def test_UD_get_corners(self):
        self.mock_station.get_xyalpha_coordinates.return_value = (.25, 3, 0)
        corners = self.detector_2.get_corners()
        self.assertEqual(corners, [(-1, 4.5), (-.5, 4.5), (-.5, 5.5), (-1, 5.5)])


class StationTests(unittest.TestCase):
    def setUp(self):
        with patch('clusters.Detector') as mock_detector:
            self.cluster = Mock()
            self.station_1 = clusters.Station(self.cluster, (0, 1), pi / 4,
                                              [(3, 4, 'LR')])
            self.mock_detector_instance = mock_detector.return_value

    def test_detector_called(self):
        with patch('clusters.Detector') as mock_detector:
            cluster = Mock()
            station = clusters.Station(cluster, (0, 1), 2, [(3, 4, 'LR')])
            mock_detector.assert_called_with(station, 3, 4, 'LR')

    def test_attributes(self):
        self.assertEqual(self.station_1.position, (0, 1))
        self.assertEqual(self.station_1.angle, pi / 4)

    def test_add_detector_for_one_instance(self):
        """
        Unfortunately, if you naively declare __detectors = [] as a *class*
        variable, you will share the same list with *all instances*.

        """
        with patch('clusters.Detector') as mock_detector:
            cluster = Mock()
            mock_detector.return_value = Mock()
            station1 = clusters.Station(cluster, (0, 1), 2, [(3, 4, 'LR')])
            mock_detector.return_value = Mock()
            station2 = clusters.Station(cluster, (0, 1), 2, [(0, 1, 'LR')])
            self.assertNotEqual(station1.detectors[0], station2.detectors[0])

    def test_detectors(self):
        self.assertEqual(self.station_1.detectors, [self.mock_detector_instance])

    def test_get_xyalpha_coordinates(self):
        with patch('clusters.Detector') as mock_detector:
            cluster = Mock()

            # Trivial
            cluster.get_xyalpha_coordinates.return_value = (0, 0, 0)
            station = clusters.Station(cluster, position=(0, 0), angle=0,
                                       detectors=[(0, 0, 'LR')])
            coordinates = station.get_xyalpha_coordinates()
            self.assertEqual(coordinates, (0, 0, 0))

            # Cluster not in origin and rotated
            cluster.get_xyalpha_coordinates.return_value = (sqrt(2) / 2, sqrt(2) / 2, pi / 8)
            station = clusters.Station(cluster, (0, 0), 0, [(0, 0, 'LR')])
            coordinates = station.get_xyalpha_coordinates()
            self.assertTupleAlmostEqual(coordinates, (sqrt(2) / 2, sqrt(2) / 2, pi / 8))

            # Station *and* cluster not in origin and cluster rotated
            cluster.get_xyalpha_coordinates.return_value = (0, 10, pi / 2)
            station = clusters.Station(cluster, (0, 5), 0, [(0, 0, 'LR')])
            coordinates = station.get_xyalpha_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10, pi / 2))

            # Station *and* cluster not in origin and cluster *and* station rotated
            cluster.get_xyalpha_coordinates.return_value = (0, 10, pi / 2)
            station = clusters.Station(cluster, (0, 5), pi / 4, [(0, 0, 'LR')])
            coordinates = station.get_xyalpha_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10, 3 * pi / 4))

    def assertTupleAlmostEqual(self, actual, expected):
        self.assertTrue(type(actual) == type(expected) == tuple)

        msg = "Tuples differ: %s != %s" % (str(actual), str(expected))
        for actual_value, expected_value in zip(actual, expected):
            self.assertAlmostEqual(actual_value, expected_value, msg=msg)


class BaseClusterTests(unittest.TestCase):
    def test_add_station(self):
        with patch('clusters.Station') as mock_station:
            cluster = clusters.BaseCluster()
            self.assertFalse(mock_station.called)

            pos = Mock(name='pos')
            angle = Mock(name='angle')
            detector_list = Mock(name='detector_list')
            cluster._add_station(pos, angle, detector_list)
            mock_station.assert_called_with(cluster, pos, angle, detector_list)

    def test_attributes(self):
        with patch('clusters.Station') as mock_station:
            mock_station_instance = Mock()
            mock_station.return_value = mock_station_instance

            cluster = clusters.BaseCluster()
            cluster._add_station(Mock(), Mock(), Mock())
            self.assertEqual(cluster.stations, [mock_station_instance])

    def test_add_station_for_one_instance(self):
        """
        Unfortunately, if you naively declare __stations = [] as a *class*
        variable, you will share the same list with *all instances*.

        """
        with patch('clusters.Station') as mock_station:
            mock_station.return_value = Mock()
            cluster1 = clusters.BaseCluster()
            cluster1._add_station((0, 0), 0, [(0, 0, 'LR')])

            mock_station.return_value = Mock()
            cluster2 = clusters.BaseCluster()
            cluster2._add_station((1, 2), 0, [(3, 0, 'LR')])

            self.assertNotEqual(cluster1.stations[0], cluster2.stations[0])

    def test_init_sets_position(self):
        cluster = clusters.BaseCluster((10., 20.), pi / 2)
        self.assertEqual(cluster._x, 10.)
        self.assertEqual(cluster._y, 20.)
        self.assertEqual(cluster._alpha, pi / 2)

    def test_get_xyalpha_coordinates(self):
        cluster = clusters.BaseCluster((10., 20.), pi / 2)
        coordinates = cluster.get_xyalpha_coordinates()
        self.assertEqual(coordinates, (10., 20., pi / 2))

    def test_set_xyalpha_coordinates(self):
        cluster = clusters.BaseCluster()
        cluster.set_xyalpha_coordinates(5., 7., pi / 2)
        self.assertEqual((cluster._x, cluster._y, cluster._alpha),
                         (5., 7., pi / 2))

    def test_set_rphialpha_coordinates(self):
        cluster = clusters.BaseCluster()
        cluster.set_rphialpha_coordinates(10., pi / 2, 0.)
        self.assertAlmostEqual(cluster._x, 0.)
        self.assertAlmostEqual(cluster._y, 10.)
        self.assertAlmostEqual(cluster._alpha, 0.)


if __name__ == '__main__':
    unittest.main()
