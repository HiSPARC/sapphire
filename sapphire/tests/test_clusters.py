from __future__ import division

from math import pi, sqrt, atan2
from numpy import array

import unittest

from mock import Mock, patch, sentinel

from sapphire import clusters


class DetectorTests(unittest.TestCase):
    def setUp(self):
        self.mock_station = Mock()
        self.detector_1 = clusters.Detector(self.mock_station, (1, 0, 0), 'LR')
        self.detector_2 = clusters.Detector(self.mock_station, (-1, 2, 1), 'UD')
        self.detector_s = clusters.Detector(self.mock_station,
                                            (sentinel.x, sentinel.y, sentinel.z),
                                            sentinel.orientation)
        self.detector_4d = clusters.Detector(station=self.mock_station,
                                             position=([0, 5], [0, 5], [0, 5]),
                                             detector_timestamps=[0, 5])

    def test_bad_arguments(self):
        self.assertRaises(Exception, clusters.Detector, self.mock_station,
                          (1, 0, 0), 'LR', [0, 1])

    def test__update_timestamp(self):
        self.assertEqual(self.detector_4d.index, -1)
        for ts, index in [(-1, 0), (0, 0), (3, 0), (7, 1), (8, 1)]:
            self.detector_4d._update_timestamp(ts)
            self.assertEqual(self.detector_4d.index, index)

    def test_4d_positions(self):
        self.mock_station.get_coordinates.return_value = (0, 0, 0, 0)
        self.assertEqual(self.detector_4d.get_coordinates(), (5, 5, 5))
        self.detector_4d._update_timestamp(3)
        self.assertEqual(self.detector_4d.get_coordinates(), (0, 0, 0))

    def test_detector_size(self):
        self.assertEqual(self.detector_1.detector_size, (.5, 1.))
        self.assertEqual(self.detector_2.detector_size, (.5, 1.))
        self.assertEqual(self.detector_s.detector_size, (.5, 1.))

    def test_get_area(self):
        self.assertEqual(self.detector_1.get_area(), .5)

    def test_attributes(self):
        self.assertIs(self.detector_1.station, self.mock_station)
        self.assertEqual(self.detector_s.x, [sentinel.x])
        self.assertEqual(self.detector_s.y, [sentinel.y])
        self.assertEqual(self.detector_s.z, [sentinel.z])
        self.assertEqual(self.detector_s.orientation, [sentinel.orientation])

    def test_get_coordinates(self):
        self.mock_station.get_coordinates.return_value = (0, 0, 0, 0)
        coordinates = self.detector_1.get_coordinates()
        self.assertEqual(coordinates, (1, 0, 0))
        xy_coordinates = self.detector_1.get_xy_coordinates()
        self.assertEqual(xy_coordinates, (1, 0))

        self.mock_station.get_coordinates.return_value = (5, 6, 0, pi / 2)
        coordinates_1 = self.detector_1.get_coordinates()
        coordinates_2 = self.detector_2.get_coordinates()
        self.assertEqual(coordinates_1, (5, 7, 0))
        self.assertEqual(coordinates_2, (3, 5, 1))
        xy_coordinates_1 = self.detector_1.get_xy_coordinates()
        xy_coordinates_2 = self.detector_2.get_xy_coordinates()
        self.assertEqual(xy_coordinates_1, (5, 7))
        self.assertEqual(xy_coordinates_2, (3, 5))

    def test_get_polar_coordinates(self):
        self.mock_station.get_coordinates.return_value = (0, 0, 0, 0)
        coordinates = self.detector_1.get_polar_coordinates()
        self.assertEqual(coordinates, (1, 0))

    def test_get_cylindrical_coordinates(self):
        self.mock_station.get_coordinates.return_value = (0, 0, 0, 0)
        coordinates = self.detector_1.get_cylindrical_coordinates()
        self.assertEqual(coordinates, (1, 0, 0))
        coordinates = self.detector_2.get_cylindrical_coordinates()
        self.assertEqual(coordinates, (sqrt(5), atan2(2, -1), 1))

    def test_LR_get_corners(self):
        self.mock_station.get_coordinates.return_value = (.25, 3, 0, 0)
        corners = self.detector_1.get_corners()
        self.assertEqual(corners, [(.75, 3.25), (.75, 2.75), (1.75, 2.75), (1.75, 3.25)])

    def test_LR_get_corners_rotated(self):
        self.mock_station.get_coordinates.return_value = (0, 0, 0, pi / 2)
        corners = self.detector_1.get_corners()
        expected_corners = [(-.25, .5), (.25, .5), (.25, 1.5), (-.25, 1.5)]
        for (x, y), (expected_x, expected_y) in zip(corners, expected_corners):
            self.assertAlmostEqual(x, expected_x)
            self.assertAlmostEqual(y, expected_y)

    def test_UD_get_corners(self):
        self.mock_station.get_coordinates.return_value = (.25, 3, 0, 0)
        corners = self.detector_2.get_corners()
        self.assertEqual(corners, [(-1, 4.5), (-.5, 4.5), (-.5, 5.5), (-1, 5.5)])

    def test_unknown_rotation_get_corners(self):
        self.mock_station.get_coordinates.return_value = (0, 0, 0, 0)
        self.detector_2.orientation = 'Foo'
        self.assertRaises(Exception, self.detector_2.get_corners)


class StationTests(unittest.TestCase):
    def setUp(self):
        with patch('sapphire.clusters.Detector') as mock_detector:
            self.cluster = Mock()
            self.station_1 = clusters.Station(self.cluster, 1, (0, 1, 2), pi / 4,
                                              [((3, 4), 'LR')])
            self.station_s = clusters.Station(self.cluster, sentinel.id,
                                              (sentinel.x, sentinel.y, sentinel.z),
                                              sentinel.angle, [],
                                              number=sentinel.number)
            self.station_4d = clusters.Station(self.cluster, 4,
                                               ([0, 5], [0, 5], [0, 5]), (0, pi),
                                               station_timestamps=[0, 5])
            self.mock_detector_instance = mock_detector.return_value

    def test_bad_arguments(self):
        with patch('sapphire.clusters.Detector'):
            self.assertRaises(Exception, clusters.Station,
                              cluster=self.cluster, station_id=1,
                              position=(0, 1, 2), station_timestamps=[1, 2])

    def test__update_timestamp(self):
        self.assertEqual(self.station_4d.index, -1)
        for ts, index in [(-1, 0), (0, 0), (3, 0), (7, 1), (8, 1)]:
            self.station_4d._update_timestamp(ts)
            self.assertEqual(self.station_4d.index, index)

    def test_4d_positions(self):
        self.cluster.get_coordinates.return_value = (0, 0, 0, 0)
        self.assertEqual(self.station_4d.get_coordinates(), (5, 5, 5, pi))
        self.station_4d._update_timestamp(3)
        self.assertEqual(self.station_4d.get_coordinates(), (0, 0, 0, 0))

    def test_get_area(self):
        self.station_1.detectors[0].get_area.return_value = .5
        self.assertEqual(self.station_1.get_area(), .5)
        for d in self.station_4d.detectors:
            d.get_area.return_value = .5
        self.assertEqual(self.station_4d.get_area([0, 1, 2]), 1.5)

    def test_detector_called(self):
        with patch('sapphire.clusters.Detector') as mock_detector:
            cluster = Mock()
            station = clusters.Station(cluster, 1, (0, 1, 2), pi, [((4, 5, 0), 'LR')])
            mock_detector.assert_called_with(station, (4, 5, 0), 'LR', [0])

    def test_attributes(self):
        self.assertEqual(self.station_s.x, [sentinel.x])
        self.assertEqual(self.station_s.y, [sentinel.y])
        self.assertEqual(self.station_s.z, [sentinel.z])
        self.assertEqual(self.station_s.angle, [sentinel.angle])
        self.assertEqual(self.station_s.number, sentinel.number)
        self.assertEqual(self.station_1.number, 1)

    def test_add_detector_for_one_instance(self):
        """
        Unfortunately, if you naively declare __detectors = [] as a *class*
        variable, you will share the same list with *all instances*.

        """
        with patch('sapphire.clusters.Detector') as mock_detector:
            cluster = Mock()
            mock_detector.return_value = Mock()
            station1 = clusters.Station(cluster, 1, (0, 1), 2, [((3, 4), 'LR')])
            mock_detector.return_value = Mock()
            station2 = clusters.Station(cluster, 2, (0, 1), 2, [((0, 1), 'LR')])
            self.assertNotEqual(station1.detectors[0], station2.detectors[0])

    def test_detectors(self):
        self.assertEqual(self.station_1.detectors, [self.mock_detector_instance])

    def test_get_coordinates(self):
        with patch('sapphire.clusters.Detector') as mock_detector:
            self.assertFalse(mock_detector.called)
            cluster = Mock()

            # Trivial
            cluster.get_coordinates.return_value = (0, 0, 0, 0)
            station = clusters.Station(cluster, 1, position=(0, 0), angle=0,
                                       detectors=[((0, 0), 'LR')])
            coordinates = station.get_coordinates()
            self.assertEqual(coordinates, (0, 0, 0, 0))
            coordinates = station.get_xyalpha_coordinates()
            self.assertEqual(coordinates, (0, 0, 0))

            # Cluster not in origin and rotated
            cluster.get_coordinates.return_value = (sqrt(2) / 2, sqrt(2) / 2, 0, pi / 8)
            station = clusters.Station(cluster, 1, (0, 0), 0, [((0, 0), 'LR')])
            coordinates = station.get_coordinates()
            self.assertTupleAlmostEqual(coordinates, (sqrt(2) / 2, sqrt(2) / 2, 0, pi / 8))
            coordinates = station.get_xyalpha_coordinates()
            self.assertTupleAlmostEqual(coordinates, (sqrt(2) / 2, sqrt(2) / 2, pi / 8))
            coordinates = station.get_xy_coordinates()
            self.assertTupleAlmostEqual(coordinates, (sqrt(2) / 2, sqrt(2) / 2))

            # Station *and* cluster not in origin and cluster rotated
            cluster.get_coordinates.return_value = (0, 10, 0, pi / 2)
            station = clusters.Station(cluster, 1, (0, 5), 0, [((0, 0), 'LR')])
            coordinates = station.get_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10, 0, pi / 2))
            coordinates = station.get_xyalpha_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10, pi / 2))
            coordinates = station.get_xy_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10))

            # Station *and* cluster not in origin and cluster *and* station rotated
            cluster.get_coordinates.return_value = (0, 10, 0, pi / 2)
            station = clusters.Station(cluster, 1, (0, 5), pi / 4, [((0, 0), 'LR')])
            coordinates = station.get_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10, 0, 3 * pi / 4))
            coordinates = station.get_xyalpha_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10, 3 * pi / 4))
            coordinates = station.get_xy_coordinates()
            self.assertTupleAlmostEqual(coordinates, (-5, 10))

            self.assertTrue(mock_detector.called)

    def test_unit_get_polar_alpha_coordinates(self):
        station_coordinates = (sqrt(2), sqrt(2), 0, pi / 4)

        cluster = Mock()
        station = clusters.Station(cluster, 1, (0, 0), 0, [])

        with patch.object(clusters.Station, 'get_coordinates') as mock_coor:
            mock_coor.return_value = station_coordinates

            r, phi, beta = station.get_polar_alpha_coordinates()
            self.assertTrue(mock_coor.called)
            self.assertAlmostEqual(r, 2)
            self.assertAlmostEqual(phi, pi / 4)
            self.assertAlmostEqual(beta, pi / 4)

    def test_get_polar_alpha_coordinates(self):
        cluster = Mock()

        # Trivial
        cluster.get_coordinates.return_value = (0, 0, 0, 0)
        station = clusters.Station(cluster, 1, position=(0, 0), angle=0,
                                   detectors=[((0, 0), 'LR')])
        coordinates = station.get_polar_alpha_coordinates()
        self.assertEqual(coordinates, (0, 0, 0))

        # Cluster not in origin and rotated
        cluster.get_coordinates.return_value = (sqrt(2) / 2, sqrt(2) / 2, 0, pi / 8)
        station = clusters.Station(cluster, 1, (0, 0), 0, [((0, 0), 'LR')])
        coordinates = station.get_polar_alpha_coordinates()
        self.assertTupleAlmostEqual(coordinates, (1, pi / 4, pi / 8))

        # Station *and* cluster not in origin and cluster rotated
        cluster.get_coordinates.return_value = (0, 10, 0, pi / 2)
        station = clusters.Station(cluster, 1, (0, 5), 0, [((0, 0), 'LR')])
        coordinates = station.get_polar_alpha_coordinates()
        self.assertTupleAlmostEqual(coordinates, (sqrt(125), 2.0344439357957027, pi / 2))

        # Station *and* cluster not in origin and cluster *and* station rotated
        cluster.get_coordinates.return_value = (0, 10, 0, pi / 2)
        station = clusters.Station(cluster, 1, (0, 5), pi / 4, [((0, 0), 'LR')])
        coordinates = station.get_polar_alpha_coordinates()
        self.assertTupleAlmostEqual(coordinates, (sqrt(125), 2.0344439357957027, 3 * pi / 4))

    def test_calc_r_and_phi_for_detectors(self):
        cluster = Mock()
        cluster.get_coordinates.return_value = (0, 0, 0, 0)
        station = clusters.Station(cluster, 1, position=(0, 0), angle=0,
                                   detectors=[((0, 0), 'LR'), ((10., 10.), 'LR')])

        r, phi = station.calc_r_and_phi_for_detectors(0, 1)
        self.assertAlmostEqual(r ** 2, 10 ** 2 + 10 ** 2)
        self.assertAlmostEqual(phi, pi / 4)

    def test_calc_center_of_mass_coordinates(self):
        cluster = Mock()
        cluster.get_coordinates.return_value = (0, 0, 0, 0)
        station = clusters.Station(cluster, 1, position=(0, 0), angle=0,
                                   detectors=[((0, 0, 0), 'LR'), ((10, 9, 1), 'LR')])
        center = station.calc_xy_center_of_mass_coordinates()
        self.assertTupleAlmostEqual(center, (5, 4.5))
        center = station.calc_center_of_mass_coordinates()
        self.assertTupleAlmostEqual(center, (5, 4.5, 0.5))

    def assertTupleAlmostEqual(self, actual, expected):
        self.assertTrue(type(actual) == type(expected) == tuple)

        msg = "Tuples differ: %s != %s" % (str(actual), str(expected))
        for actual_value, expected_value in zip(actual, expected):
            self.assertAlmostEqual(actual_value, expected_value, msg=msg)


class BaseClusterTests(unittest.TestCase):
    def test_add_station(self):
        with patch('sapphire.clusters.Station') as mock_station:
            cluster = clusters.BaseCluster()
            self.assertFalse(mock_station.called)

            x = Mock(name='x')
            y = Mock(name='y')
            z = Mock(name='z')
            angle = Mock(name='angle')
            detector_list = Mock(name='detector_list')
            number = Mock(name='number')
            cluster._add_station((x, y, z), angle, detector_list, number=number)
            mock_station.assert_called_with(cluster, 0, (x, y, z), angle,
                                            detector_list, [0], [0], number)

    def test_set_timestamp(self):
        with patch('sapphire.clusters.Station'):
            cluster = clusters.BaseCluster()
            cluster._add_station(sentinel.position)
            self.assertEqual(cluster._timestamp, 2147483647)
            cluster.set_timestamp(sentinel.timestamp)
            self.assertEqual(cluster._timestamp, sentinel.timestamp)
            cluster.stations[0]._update_timestamp.assert_called_with(sentinel.timestamp)

    def test_attributes(self):
        with patch('sapphire.clusters.Station') as mock_station:
            mock_station_instance = Mock()
            mock_station.return_value = mock_station_instance

            cluster = clusters.BaseCluster()
            cluster._add_station(Mock(), Mock(), Mock(), Mock())
            self.assertEqual(cluster.stations, [mock_station_instance])

    def test_add_station_for_one_instance(self):
        """
        Unfortunately, if you naively declare __stations = [] as a *class*
        variable, you will share the same list with *all instances*.

        """
        with patch('sapphire.clusters.Station') as mock_station:
            mock_station.return_value = Mock()
            cluster1 = clusters.BaseCluster()
            cluster1._add_station((0, 0, 0), 0, [(0, 0, 'LR')])

            mock_station.return_value = Mock()
            cluster2 = clusters.BaseCluster()
            cluster2._add_station((1, 2, 0), 0, [(3, 0, 'LR')])

            self.assertNotEqual(cluster1.stations[0], cluster2.stations[0])

    def test_get_station_by_number(self):
        cluster = clusters.BaseCluster((0, 0, 0), 0)
        cluster._add_station((0, 0, 0), 0, number=501)
        self.assertEqual(cluster.get_station(501), cluster.stations[0])

    def test_init_sets_position(self):
        cluster = clusters.BaseCluster((10., 20.), pi / 2)
        self.assertEqual(cluster.x, 10.)
        self.assertEqual(cluster.y, 20.)
        self.assertEqual(cluster.z, 0.)
        self.assertEqual(cluster.alpha, pi / 2)

    def test_get_coordinates(self):
        cluster = clusters.BaseCluster((10., 20., 0), pi / 2)
        coordinates = cluster.get_coordinates()
        self.assertEqual(coordinates, (10., 20., 0, pi / 2))
        coordinates = cluster.get_xyalpha_coordinates()
        self.assertEqual(coordinates, (10., 20., pi / 2))
        coordinates = cluster.get_xy_coordinates()
        self.assertEqual(coordinates, (10., 20.))

    def test_get_polar_alpha_coordinates(self):
        cluster = clusters.BaseCluster((-sqrt(2) / 2, sqrt(2) / 2), pi / 2)
        r, phi, alpha = cluster.get_polar_alpha_coordinates()
        self.assertAlmostEqual(r, 1.)
        self.assertAlmostEqual(phi, 3 * pi / 4)
        self.assertEqual(alpha, pi / 2)

    def test_set_coordinates(self):
        cluster = clusters.BaseCluster()
        cluster.set_coordinates(sentinel.x, sentinel.y, sentinel.z, sentinel.alpha)
        self.assertEqual((cluster.x, cluster.y, cluster.z, cluster.alpha),
                         (sentinel.x, sentinel.y, sentinel.z, sentinel.alpha))

    def test_set_rphialpha_coordinates(self):
        cluster = clusters.BaseCluster()
        cluster.set_cylindrical_coordinates(10., pi / 2, sentinel.z, sentinel.alpha)
        self.assertAlmostEqual(cluster.x, 0.)
        self.assertAlmostEqual(cluster.y, 10.)
        self.assertAlmostEqual(cluster.z, sentinel.z)
        self.assertAlmostEqual(cluster.alpha, sentinel.alpha)

    def test_calc_r_and_phi_for_stations(self):
        cluster = clusters.BaseCluster()
        cluster._add_station((0, 0), 0)
        cluster._add_station((1, sqrt(3)), 0)
        r, phi, z = cluster.calc_rphiz_for_stations(0, 1)
        self.assertAlmostEqual(r, 2)
        self.assertAlmostEqual(phi, pi / 3.)
        self.assertAlmostEqual(z, 0)

    def test_calc_xy_center_of_mass_coordinates(self):
        cluster = clusters.BaseCluster()
        cluster._add_station((0, 0), 0, [((0, 5 * sqrt(3)), 'UD'),
                                         ((0, 5 * sqrt(3) / 3), 'UD'),
                                         ((-10, 0), 'LR'),
                                         ((10, 0), 'LR')])

        x, y = cluster.calc_xy_center_of_mass_coordinates()
        self.assertAlmostEqual(x, 0)
        self.assertAlmostEqual(y, 5 * sqrt(3) / 3)

    def test__distance(self):
        x = array([-5., 4., 3.])
        y = array([2., -1., 0.])
        dist = clusters.BaseCluster()._distance(x, y)
        self.assertAlmostEqual(dist, sqrt(49 + 25 + 9))
        x = array([-5., 4.])
        y = array([2., -1.])
        dist = clusters.BaseCluster()._distance(x, y)
        self.assertAlmostEqual(dist, sqrt(49 + 25))

    def test_calc_distance_between_stations(self):
        cluster = clusters.BaseCluster()
        cluster._add_station((0, 0, 0), 0, number=1)
        cluster._add_station((3, 4, 5), 0, number=2)
        dist = cluster.calc_distance_between_stations(1, 2)
        self.assertAlmostEqual(dist, sqrt(9 + 16 + 25))
        dist = cluster.calc_distance_between_stations(1, 0)
        self.assertIsNone(dist)

    def test_calc_horizontal_distance_between_stations(self):
        cluster = clusters.BaseCluster()
        cluster._add_station((0, 0, 0), 0, number=1)
        cluster._add_station((3, 4, 5), 0, number=2)
        dist = cluster.calc_horizontal_distance_between_stations(1, 2)
        self.assertAlmostEqual(dist, sqrt(9 + 16))
        dist = cluster.calc_distance_between_stations(1, 0)
        self.assertIsNone(dist)


class CompassStationsTests(unittest.TestCase):
    def test_cluster_stations(self):
        cluster = clusters.CompassStations()
        cluster._add_station((0, 0, 0), [(5, 270, 0, 0)])
        stations = cluster.stations
        self.assertEqual(len(stations), 1)
        detectors = stations[0].detectors
        self.assertEqual(len(detectors), 1)
        self.assertEqual(detectors[0].x, [-5.0])
        self.assertTupleAlmostEqual(detectors[0].get_coordinates(),
                                    (-5.0, 0, 0))

    def assertTupleAlmostEqual(self, actual, expected):
        self.assertTrue(type(actual) == type(expected) == tuple)

        msg = "Tuples differ: %s != %s" % (str(actual), str(expected))
        for actual_value, expected_value in zip(actual, expected):
            self.assertAlmostEqual(actual_value, expected_value, msg=msg)


class SimpleClusterTests(unittest.TestCase):
    def test_init_calls_super_init(self):
        with patch.object(clusters.BaseCluster, '__init__',
                          mocksignature=True) as mock_base_init:
            clusters.SimpleCluster()
            self.assertTrue(mock_base_init.called)

    def test_get_coordinates_after_init(self):
        cluster = clusters.SimpleCluster()
        coordinates = cluster.get_xyalpha_coordinates()
        self.assertEqual(coordinates, (0., 0., 0.))

    def test_cluster_stations(self):
        cluster = clusters.SimpleCluster()
        stations = cluster.stations
        self.assertEqual(len(stations), 4)


class SingleStationTests(unittest.TestCase):
    def test_init_calls_super_init(self):
        with patch.object(clusters.BaseCluster, '__init__',
                          mocksignature=True) as mock_base_init:
            clusters.SingleStation()
            self.assertTrue(mock_base_init.called)

    def test_get_coordinates_after_init(self):
        cluster = clusters.SingleStation()
        coordinates = cluster.get_xyalpha_coordinates()
        self.assertEqual(coordinates, (0., 0., 0.))

    def test_single_station(self):
        cluster = clusters.SingleStation()
        stations = cluster.stations
        self.assertEqual(len(stations), 1)


class SingleDetectorStationTests(unittest.TestCase):
    def setUp(self):
        self.cluster = clusters.SingleDetectorStation()

    def test_init_calls_super_init(self):
        with patch.object(clusters.BaseCluster, '__init__',
                          mocksignature=True) as mock_base_init:
            clusters.SingleDetectorStation()
            self.assertTrue(mock_base_init.called)

    def test_get_coordinates_after_init(self):
        coordinates = self.cluster.get_xyalpha_coordinates()
        self.assertEqual(coordinates, (0., 0., 0.))

    def test_single_station_single_detector(self):
        self.assertEqual(len(self.cluster.stations), 1)
        self.assertEqual(len(self.cluster.stations[0].detectors), 1)


class SingleTwoDetectorStationTests(unittest.TestCase):
    def test_init_calls_super_init(self):
        with patch.object(clusters.BaseCluster, '__init__',
                          mocksignature=True) as mock_base_init:
            clusters.SingleTwoDetectorStation()
            self.assertTrue(mock_base_init.called)

    def test_get_coordinates_after_init(self):
        cluster = clusters.SingleTwoDetectorStation()
        coordinates = cluster.get_xyalpha_coordinates()
        self.assertEqual(coordinates, (0., 0., 0.))

    def test_single_station(self):
        cluster = clusters.SingleTwoDetectorStation()
        stations = cluster.stations
        self.assertEqual(len(stations), 1)
        detectors = stations[0].detectors
        self.assertEqual(len(detectors), 2)


class SingleDiamondStationTests(unittest.TestCase):
    def test_init_calls_super_init(self):
        with patch.object(clusters.BaseCluster, '__init__',
                          mocksignature=True) as mock_base_init:
            clusters.SingleDiamondStation()
            self.assertTrue(mock_base_init.called)

    def test_get_coordinates_after_init(self):
        cluster = clusters.SingleDiamondStation()
        coordinates = cluster.get_xyalpha_coordinates()
        self.assertEqual(coordinates, (0., 0., 0.))

    def test_single_station(self):
        cluster = clusters.SingleDiamondStation()
        stations = cluster.stations
        self.assertEqual(len(stations), 1)
        detectors = stations[0].detectors
        self.assertEqual(len(detectors), 4)


if __name__ == '__main__':
    unittest.main()
