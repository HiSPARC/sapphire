import unittest
import warnings

from mock import sentinel, patch, Mock, MagicMock
from numpy import isnan, nan, pi, sqrt, arcsin, arctan

from sapphire.analysis import direction_reconstruction


class EventDirectionReconstructionTest(unittest.TestCase):

    def test_init(self):
        dirrec = direction_reconstruction.EventDirectionReconstruction(sentinel.station)
        self.assertEqual(dirrec.direct, direction_reconstruction.DirectAlgorithmCartesian3D)
        self.assertEqual(dirrec.fit, direction_reconstruction.RegressionAlgorithm3D)
        self.assertEqual(dirrec.station, sentinel.station)

    def test_set_cluster_timestamp(self):
        station = Mock()
        dirrec = direction_reconstruction.EventDirectionReconstruction(station)
        theta, phi, ids = dirrec.reconstruct_event({'timestamp': sentinel.timestamp}, detector_ids=[])
        station.cluster.set_timestamp.assert_called_with(sentinel.timestamp)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))

    @patch.object(direction_reconstruction, 'detector_arrival_time')
    def test_bad_times(self, mock_detector_arrival_time):
        mock_detector_arrival_time.return_value = nan
        station = Mock()
        dirrec = direction_reconstruction.EventDirectionReconstruction(station)
        theta, phi, ids = dirrec.reconstruct_event({'timestamp': sentinel.timestamp}, detector_ids=[0, 1])
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))

    @patch.object(direction_reconstruction, 'detector_arrival_time')
    def test_reconstruct_event(self, mock_detector_arrival_time):
        mock_detector_arrival_time.return_value = 0.
        station = MagicMock()
        detector = Mock()
        detector.get_coordinates.return_value = [sentinel.x, sentinel.y, sentinel.z]
        station.detectors.__getitem__.side_effect = lambda name: detector
        dirrec = direction_reconstruction.EventDirectionReconstruction(station)
        dirrec.direct = Mock()
        dirrec.fit = Mock()
        dirrec.direct.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        dirrec.fit.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        event = {'timestamp': sentinel.timestamp}

        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=[0, 1])
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 0)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(ids), 2)

        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=[0, 1, 2])
        dirrec.direct.reconstruct_common.assert_called_once_with(
            [0.] * 3, [sentinel.x] * 3, [sentinel.y] * 3, [sentinel.z] * 3)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 0)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(ids), 3)

        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=[0, 1, 2, 3])
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        dirrec.fit.reconstruct_common.assert_called_once_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(ids), 4)
        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=None)
        dirrec.fit.reconstruct_common.assert_called_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 2)

    @patch.object(direction_reconstruction.EventDirectionReconstruction, 'reconstruct_event')
    def test_reconstruct_events(self, mock_reconstruct_event):
        mock_reconstruct_event.return_value = [sentinel.theta, sentinel.phi, sentinel.ids]
        dirrec = direction_reconstruction.EventDirectionReconstruction(sentinel.station)
        self.assertEqual(dirrec.reconstruct_events([sentinel.event, sentinel.event], sentinel.detector_ids, sentinel.offsets, progress=False),
                         ((sentinel.theta, sentinel.theta), (sentinel.phi, sentinel.phi), (sentinel.ids, sentinel.ids)))
        self.assertEqual(mock_reconstruct_event.call_count, 2)
        mock_reconstruct_event.assert_called_with(sentinel.event, sentinel.detector_ids, sentinel.offsets)
        self.assertEqual(dirrec.reconstruct_events([], sentinel.detector_ids, sentinel.offsets, progress=False),
                         ((), (), ()))
        self.assertEqual(mock_reconstruct_event.call_count, 2)


class CoincidenceDirectionReconstructionTest(unittest.TestCase):

    def setUp(self):
        self.dirrec = direction_reconstruction.CoincidenceDirectionReconstruction(sentinel.cluster)

    def test_init(self):
        dirrec = direction_reconstruction.CoincidenceDirectionReconstruction(sentinel.cluster)
        self.assertEqual(dirrec.direct, direction_reconstruction.DirectAlgorithmCartesian3D)
        self.assertEqual(dirrec.fit, direction_reconstruction.RegressionAlgorithm3D)
        self.assertEqual(dirrec.cluster, sentinel.cluster)

    def test_set_cluster_timestamp(self):
        dirrec = self.dirrec
        cluster = Mock()
        dirrec.cluster = cluster
        coincidence = [[sentinel.station_number, {'timestamp': 1}], [0, 0], [0, 0]]
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence, station_numbers=[])
        cluster.set_timestamp.assert_called_with(1)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))

    @patch.object(direction_reconstruction, 'station_arrival_time')
    def test_reconstruct_coincidence(self, mock_station_arrival_time):
        dirrec = self.dirrec
        mock_station_arrival_time.return_value = 0.
        cluster = MagicMock()
        station = MagicMock()
        cluster.get_station.return_value = station
        station.calc_center_of_mass_coordinates.return_value = [sentinel.x, sentinel.y, sentinel.z]
        dirrec.cluster = cluster
        dirrec.direct = Mock()
        dirrec.fit = Mock()
        dirrec.direct.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        dirrec.fit.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        coincidence_2 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event]]
        coincidence_3 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event]]
        coincidence_4 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event], [3, sentinel.event]]

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_2)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 0)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_3, station_numbers=[1, 2])
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 0)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 2)

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_3)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.direct.reconstruct_common.assert_called_once_with(
            [0.] * 3, [sentinel.x] * 3, [sentinel.y] * 3, [sentinel.z] * 3)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 0)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 3)

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.fit.reconstruct_common.assert_called_once_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 4)

        mock_station_arrival_time.return_value = nan

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4)
        cluster.set_timestamp.assert_called_with(1)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 1)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

    @patch.object(direction_reconstruction.CoincidenceDirectionReconstruction, 'reconstruct_coincidence')
    def test_reconstruct_coincidences(self, mock_reconstruct_coincidence):
        mock_reconstruct_coincidence.return_value = [sentinel.theta, sentinel.phi, sentinel.nums]
        self.assertEqual(self.dirrec.reconstruct_coincidences([sentinel.coincidence, sentinel.coincidence], sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((sentinel.theta, sentinel.theta), (sentinel.phi, sentinel.phi), (sentinel.nums, sentinel.nums)))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)
        mock_reconstruct_coincidence.assert_called_with(sentinel.coincidence, sentinel.station_numbers, sentinel.offsets)
        self.assertEqual(self.dirrec.reconstruct_coincidences([], sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((), (), ()))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)


class CoincidenceDirectionReconstructionDetectorsTest(CoincidenceDirectionReconstructionTest):

    def setUp(self):
        self.dirrec = direction_reconstruction.CoincidenceDirectionReconstructionDetectors(sentinel.cluster)

    @patch.object(direction_reconstruction, 'relative_detector_arrival_times')
    def test_reconstruct_coincidence(self, mock_arrival_times):
        dirrec = self.dirrec
        mock_arrival_times.return_value = [0., nan, nan, nan]
        cluster = MagicMock()
        station = MagicMock()
        cluster.get_station.return_value = station
        detector = MagicMock()
        detector.get_coordinates.return_value = [sentinel.x, sentinel.y, sentinel.z]
        station.detectors = [detector] * 4
        dirrec.cluster = cluster
        dirrec.direct = Mock()
        dirrec.fit = Mock()
        dirrec.direct.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        dirrec.fit.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        coincidence_2 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event]]
        coincidence_3 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event]]
        coincidence_4 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event], [3, sentinel.event]]

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_2)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 0)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_3, station_numbers=[1, 2])
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 0)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 0)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 2)

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_3)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.direct.reconstruct_common.assert_called_once_with(
            [0.] * 3, [sentinel.x] * 3, [sentinel.y] * 3, [sentinel.z] * 3)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 0)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 3)

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.fit.reconstruct_common.assert_called_once_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 4)

        mock_arrival_times.return_value = [nan] * 4

        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4)
        cluster.set_timestamp.assert_called_with(1)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 1)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

    @patch.object(direction_reconstruction.CoincidenceDirectionReconstructionDetectors, 'reconstruct_coincidence')
    def test_reconstruct_coincidences(self, mock_reconstruct_coincidence):
        mock_reconstruct_coincidence.return_value = [sentinel.theta, sentinel.phi, sentinel.nums]
        self.assertEqual(self.dirrec.reconstruct_coincidences([sentinel.coincidence, sentinel.coincidence], sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((sentinel.theta, sentinel.theta), (sentinel.phi, sentinel.phi), (sentinel.nums, sentinel.nums)))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)
        mock_reconstruct_coincidence.assert_called_with(sentinel.coincidence, sentinel.station_numbers, sentinel.offsets)
        self.assertEqual(self.dirrec.reconstruct_coincidences([], sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((), (), ()))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)


class BaseAlgorithm(object):

    """Use this class to check the different algorithms

    They should give similar results and errors in some cases.
    Theses tests use three detections at same height.

    """

    def call_reconstruct(self, t, x, y, z):
        return self.algorithm.reconstruct_common(t, x, y, z)

    def test_stations_in_line(self):
        """Three detection points on a line do not provide a solution."""

        # On a line in x
        t = (0., 2., 3.)
        x = (0., 0., 0.)  # same x
        y = (0., 5., 10.)
        z = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

        # Diagonal line
        t = (0., 2., 3.)
        x = (0., 5., 10.)
        y = (0., 5., 10.)
        z = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

    def test_same_stations(self):
        """Multiple detections at same point make reconstruction impossible.

        With different arrival time.

        """
        # Two at same location
        t = (0., 2., 3.)
        x = (0., 0., 1.)
        y = (5., 5., 6.)
        z = (0., 0., 1.)
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

        t = (0., 2., 3.)
        x = (0., 1., 0.)
        y = (5., 6., 5.)
        z = (0., 1., 0.)
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

        t = (0., 2., 3.)
        x = (1., 0., 0.)
        y = (6., 5., 5.)
        z = (1., 0., 0.)
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

        # Three at same location
        t = (0., 2., 3.)
        x = (0., 0., 0.)  # same x
        y = (5., 5., 5.)  # same y
        z = (0., 0., 0.)  # same z
        result = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(result).all())

    def test_shower_from_above(self):
        """Simple shower from zenith, azimuth can be any allowed value."""

        t = (0., 0., 0.)  # same t
        x = (0., 10., 0.)
        y = (0., 0., 10.)
        z = (0., 0., 0.)  # same z
        theta, phi = self.call_reconstruct(t, x, y, z)
        self.assertAlmostEqual(theta, 0., 4)
        # azimuth can be any value between -pi and pi
        self.assertTrue(-pi <= phi <= pi)

    def test_to_large_dt(self):
        """Time difference larger than expected by speed of light."""

        # TODO: Add better test with smaller tolerance

        x = (0., -5., 5.)
        y = (sqrt(100 - 25), 0., 0.)
        z = (0., 0., 0.)

        t = (35., 0., 0.)
        theta, phi = self.call_reconstruct(t, x, y, z)
        self.assertTrue(isnan(theta))

    def test_showers_at_various_angles(self):
        """Simple shower from specific zenith angles."""

        c = 0.299792458

        x = (0., -5., 5.)
        y = (sqrt(100 - 25), 0., 0.)
        z = (0., 0., 0.)

        # triangle height
        h = sqrt(100 - 25)

        times = (2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5)

        for time in times:
            for i in range(3):
                zenith = arcsin((time * c) / h)

                t = [0., 0., 0.]
                t[i] = time
                azimuths = [-pi / 2, pi / 6, pi * 5 / 6]
                theta, phi = self.call_reconstruct(t, x, y, z)
                self.assertAlmostEqual(phi, azimuths[i], 4)
                self.assertAlmostEqual(theta, zenith, 4)
                # Compare with z=None, should default to 0
                theta_no_z, phi_no_z = self.call_reconstruct(t, x, y, None)
                self.assertEqual((theta, phi), (theta_no_z, phi_no_z))

                t = [time] * 3
                t[i] = 0.
                azimuths = [pi / 2, -pi * 5 / 6, -pi / 6]
                theta, phi = self.call_reconstruct(t, x, y, z)
                self.assertAlmostEqual(phi, azimuths[i], 4)
                self.assertAlmostEqual(theta, zenith, 4)
                # Compare with z=None, should default to 0
                theta_no_z, phi_no_z = self.call_reconstruct(t, x, y, None)
                self.assertEqual((theta, phi), (theta_no_z, phi_no_z))


class DirectAlgorithm(BaseAlgorithm):

    """Use this class to check algorithms that only support three detections

    They should give similar warnings in some cases.

    """

    def test_to_many_stations(self):
        """To many stations should issue a warning.

        Moreover, the result should be based on the first three detections

        """
        # Shower from above (for first three detectors)
        x = (0., 10., 0., 10.)
        y = (0., 0., 10., 10.)
        z = (0., 0., 0., 0.)
        t = (0., 0., 0., 10.)

        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            theta, phi = self.call_reconstruct(t, x, y, z)
            self.assertTrue(issubclass(w[0].category, UserWarning))
            self.assertAlmostEqual(theta, 0., 4)
            self.assertTrue(-pi <= phi <= pi)


class AltitudeAlgorithm(BaseAlgorithm):

    """Use this class to check the altitude support

    They should give similar results and errors in some cases.

    """

    def test_stations_altitude(self):
        """Simple shower on a non horizontal square."""

        x = (0., 10., 10.)
        y = (0, 0., 10.)
        z = (2., 0., -2.)

        zenith = arctan(4. / 10. / sqrt(2))

        t = [0., 0., 0.]
        azimuth = pi / 4.
        theta, phi = self.call_reconstruct(t, x, y, z)

        self.assertAlmostEqual(phi, azimuth, 5)
        self.assertAlmostEqual(theta, zenith, 5)


class DirectAltitudeAlgorithm(DirectAlgorithm, AltitudeAlgorithm):

    """Test algorithm that uses only 3 detectors and has altitude support."""

    pass


class MultiAlgorithm(BaseAlgorithm):

    """Use this class to check the different algorithms for more stations

    They should give similar results and errors in some cases.

    """

    def test_diamond_stations(self):
        """Simple shower from specific zenith angles."""

        c = 0.299792458

        x = (0., -5., 5., 10.)
        y = (sqrt(100 - 25), 0., 0., sqrt(100 - 25))
        z = (0., 0., 0., 0.)

        # triangle height
        h = sqrt(100 - 25)

        times = (2.5, 5., 7.5, 10., 12.5, 15., 17.5, 20., 22.5, 25., 27.5)

        for time in times:
            zenith = arcsin((time * c) / h)

            t = [0., 0., 0., 0.]
            t[1] = time
            t[3] = -time
            azimuth = pi / 6
            theta, phi = self.call_reconstruct(t, x, y, z)
            self.assertAlmostEqual(phi, azimuth, 5)
            self.assertAlmostEqual(theta, zenith, 5)

    def test_square_stations(self):
        """Simple shower from specific zenith angles."""

        c = 0.299792458

        x = (0., 5., 5., 0.)
        y = (0, 0., 5., 5.)
        z = (0., 0., 0., 0.)

        # triangle height
        h = sqrt(50. / 4.)

        times = (2.5, 5., 7.5, 10.)

        for time in times:
            zenith = arcsin((time * c) / h)

            t = [0., 0., 0., 0.]
            t[0] = -time
            t[2] = time
            azimuth = - 3 * pi / 4
            theta, phi = self.call_reconstruct(t, x, y, z)
            self.assertAlmostEqual(phi, azimuth, 5)
            self.assertAlmostEqual(theta, zenith, 5)


class MultiAltitudeAlgorithm(MultiAlgorithm, AltitudeAlgorithm):

    """Check some algorithms for multiple stations at different altitudes.

    They should give similar results and errors in some cases.

    """

    def test_hexagon_altitude(self):
        """Simple shower on a non horizontal square."""

        x = (-5., 5., 10., 5., -5., -10.)
        y = (-5. * sqrt(3), -5. * sqrt(3), 0., 5. * sqrt(3), 5. * sqrt(3), 0.)
        z = (0., -3., -5., -3., 0., 4.)

        zenith = 0.38333
        azimuth = 0.00000

        t = [0., 0., 0., 0., 0., 0.]
        theta, phi = self.call_reconstruct(t, x, y, z)

        self.assertAlmostEqual(phi, azimuth, 4)
        self.assertAlmostEqual(theta, zenith, 4)


class DirectAlgorithmTest(unittest.TestCase, DirectAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithm()


class DirectAlgorithmCartesian2DTest(unittest.TestCase, DirectAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian2D()


class DirectAlgorithmCartesian3DTest(unittest.TestCase,
                                     DirectAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian3D()


class FitAlgorithmTest(unittest.TestCase, MultiAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.FitAlgorithm()


class RegressionAlgorithmTest(unittest.TestCase, MultiAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.RegressionAlgorithm()


class RegressionAlgorithm3DTest(unittest.TestCase, MultiAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.RegressionAlgorithm3D()


if __name__ == '__main__':
    unittest.main()
