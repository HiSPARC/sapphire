import unittest
import warnings

from mock import sentinel, patch, Mock, MagicMock
from numpy import isnan, nan, pi, sqrt, arcsin, arctan, array

from sapphire.analysis import direction_reconstruction
from sapphire.simulations.showerfront import ConeFront


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

        # To few detectors
        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=[0, 1])
        dirrec.direct.reconstruct_common.assert_not_called()
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(ids), 2)

        # Three detections, direct reconstruction
        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=[0, 1, 2])
        dirrec.direct.reconstruct_common.assert_called_once_with(
            [0.] * 3, [sentinel.x] * 3, [sentinel.y] * 3, [sentinel.z] * 3, None)
        dirrec.fit.reconstruct_common.assert_not_called()
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(ids), 3)

        # Four detections, fit reconstruction
        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=[0, 1, 2, 3])
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        dirrec.fit.reconstruct_common.assert_called_once_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4, None)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(ids), 4)
        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=None)
        dirrec.fit.reconstruct_common.assert_called_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4, None)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 2)

        # Four detections, fit reconstruction with offsets
        offsets = MagicMock(spec=direction_reconstruction.Station)
        offsets.detector_timing_offset.return_value = [sentinel.offset] * 4
        theta, phi, ids = dirrec.reconstruct_event(event, detector_ids=[0, 1, 2, 3], offsets=offsets)
        offsets.detector_timing_offset.assert_called_once_with(sentinel.timestamp)
        mock_detector_arrival_time.assert_called_with(event, 3, offsets.detector_timing_offset.return_value)

    @patch.object(direction_reconstruction.EventDirectionReconstruction, 'reconstruct_event')
    def test_reconstruct_events(self, mock_reconstruct_event):
        mock_reconstruct_event.return_value = [sentinel.theta, sentinel.phi, sentinel.ids]
        dirrec = direction_reconstruction.EventDirectionReconstruction(sentinel.station)
        self.assertEqual(dirrec.reconstruct_events([sentinel.event, sentinel.event],
                                                   sentinel.detector_ids, sentinel.offsets, progress=False),
                         ((sentinel.theta, sentinel.theta), (sentinel.phi, sentinel.phi), (sentinel.ids, sentinel.ids)))
        self.assertEqual(mock_reconstruct_event.call_count, 2)
        mock_reconstruct_event.assert_called_with(sentinel.event, sentinel.detector_ids, sentinel.offsets, None)
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
        dirrec.curved = Mock()
        dirrec.direct.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        dirrec.fit.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        dirrec.curved.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        coincidence_2 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event]]
        coincidence_3 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event]]
        coincidence_4 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event], [3, sentinel.event]]

        # To few events
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_2)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

        # To few eligible events
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_3, station_numbers=[1, 2])
        dirrec.direct.reconstruct_common.assert_not_called()
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 2)

        # Three events, no initial core, direct reconstruction
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_3)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.direct.reconstruct_common.assert_called_once_with(
            [0.] * 3, [sentinel.x] * 3, [sentinel.y] * 3, [sentinel.z] * 3, {})
        dirrec.fit.reconstruct_common.assert_not_called()
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 3)

        # Four events, no initial core, fit reconstruction
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.fit.reconstruct_common.assert_called_once_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4, {})
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        dirrec.curved.reconstruct_common.assert_not_called()
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 4)

        # Four events, with initial core, curved reconstruction
        initial = {'core_x': sentinel.core_x, 'core_y': sentinel.core_y}
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4, initial=initial)
        cluster.set_timestamp.assert_called_with(1)
        self.assertEqual(dirrec.curved.reconstruct_common.call_count, 1)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 4)

        # Four events, no valid station arrival times
        mock_station_arrival_time.return_value = nan
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4)
        cluster.set_timestamp.assert_called_with(1)
        self.assertEqual(dirrec.curved.reconstruct_common.call_count, 1)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 1)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

    @patch.object(direction_reconstruction.CoincidenceDirectionReconstruction, 'reconstruct_coincidence')
    def test_reconstruct_coincidences(self, mock_reconstruct_coincidence):
        mock_reconstruct_coincidence.return_value = [sentinel.theta, sentinel.phi, sentinel.nums]
        self.assertEqual(self.dirrec.reconstruct_coincidences([sentinel.coincidence, sentinel.coincidence],
                                                              sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((sentinel.theta, sentinel.theta), (sentinel.phi, sentinel.phi), (sentinel.nums, sentinel.nums)))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)
        mock_reconstruct_coincidence.assert_called_with(sentinel.coincidence, sentinel.station_numbers, sentinel.offsets, None)
        self.assertEqual(self.dirrec.reconstruct_coincidences([], sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((), (), ()))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)

    def test_get_station_offsets(self):
        dirrec = self.dirrec
        mock_offsets = Mock()
        mock_offsets.return_value = sentinel.best_offset
        dirrec.determine_best_offsets = mock_offsets
        coincidence_events = [(sentinel.sn1, None)]
        station_numbers = None
        offsets = {}
        ts0 = 86400
        result = dirrec.get_station_offsets(coincidence_events, station_numbers,
                                            offsets, ts0)
        self.assertEqual(result, offsets)

        offsets = {1: MagicMock(spec=direction_reconstruction.Station)}
        result = dirrec.get_station_offsets(coincidence_events, station_numbers,
                                            offsets, ts0)
        self.assertEqual(result, sentinel.best_offset)
        mock_offsets.assert_called_once_with([sentinel.sn1], ts0, offsets)

        ts0 = 864000 + 12345
        result = dirrec.get_station_offsets(coincidence_events, station_numbers,
                                            offsets, ts0)
        self.assertEqual(result, sentinel.best_offset)
        mock_offsets.assert_called_with([sentinel.sn1], 864000, offsets)

        station_numbers = sentinel.station_numbers
        result = dirrec.get_station_offsets(coincidence_events, station_numbers,
                                            offsets, ts0)
        self.assertEqual(result, sentinel.best_offset)
        mock_offsets.assert_called_with(sentinel.station_numbers, 864000, offsets)

    def test_determine_best_offsets(self):
        dirrec = self.dirrec
        mock_offsets = Mock()
        mock_offsets.station_timing_offset.return_value = (1, 2)
        mock_offsets.detector_timing_offset.return_value = [1, 0, 2, 3]
        offsets = {sn: mock_offsets for sn in [1, 2, 3]}
        station_numbers = [1, 2]
        midnight_ts = sentinel.midnight_ts
        best_offsets = dirrec.determine_best_offsets(station_numbers, midnight_ts, offsets)
        self.assertEqual(list(best_offsets.keys()), station_numbers)
        self.assertEqual(list(best_offsets.values()), [[1.0, 0.0, 2.0, 3.0],
                                                       [2.0, 1.0, 3.0, 4.0]])

    def test_determine_best_reference(self):
        # last station would be best reference, but not in station_numbers
        # second and third station are tied, so second is best reference
        error_matrix = array([[0, 5, 2, 1],
                              [5, 0, 1, 1],
                              [2, 1, 0, 1],
                              [1, 1, 1, 0]])
        station_numbers = [1, 2, 3]
        ref, pred = self.dirrec.determine_best_reference(error_matrix, station_numbers)
        self.assertEqual(ref, 2)
        predecessors = array([[-9999, 3, 0, 0],
                              [3, -9999, 1, 1],
                              [2, 2, -9999, 2],
                              [3, 3, 3, -9999]])
        self.assertEqual(pred.tolist(), predecessors.tolist())

    def test__reconstruct_best_offset(self):
        offset = self.dirrec._reconstruct_best_offset([], 1, 1, [], [])
        self.assertEqual(offset, 0)

        predecessors = array([[-9999, 0, 1],
                              [1, -9999, 1],
                              [1, 2, -9999]])
        offset_matrix = array([[0, -1, -1],
                               [1, 0, -1],
                               [1, 1, 0]])
        station_numbers = [1, 2, 3]

        combinations = [(1, 1, 0),
                        (1, 2, 1),
                        (1, 3, 2),
                        (2, 3, 1)]
        for sn1, sn2, offset in combinations:
            o12 = self.dirrec._reconstruct_best_offset(predecessors, sn1, sn2, station_numbers, offset_matrix)
            o21 = self.dirrec._reconstruct_best_offset(predecessors, sn2, sn1, station_numbers, offset_matrix)
            self.assertEqual(offset, o12)
            self.assertEqual(o21, -o12)

    def test__calculate_offsets(self):
        mock_station = Mock()
        mock_station.detector_timing_offset.return_value = [0., 1., 2., 3.]
        offset = 1.
        ts0 = sentinel.timestamp
        offsets = self.dirrec._calculate_offsets(mock_station, ts0, offset)
        mock_station.detector_timing_offset.assert_called_once_with(ts0)
        self.assertEqual(offsets, [1., 2., 3., 4.])


class CoincidenceDirectionReconstructionDetectorsTest(CoincidenceDirectionReconstructionTest):

    def setUp(self):
        self.dirrec = direction_reconstruction.CoincidenceDirectionReconstructionDetectors(sentinel.cluster)

    @patch.object(direction_reconstruction, 'relative_detector_arrival_times')
    def test_reconstruct_coincidence(self, mock_arrival_times):
        dirrec = self.dirrec
        mock_arrival_times.return_value = [0., 0., nan, nan]
        cluster = MagicMock()
        station = MagicMock()
        cluster.get_station.return_value = station
        detector = MagicMock()
        detector.get_coordinates.return_value = [sentinel.x, sentinel.y, sentinel.z]
        station.detectors = [detector] * 4
        dirrec.cluster = cluster
        dirrec.direct = Mock()
        dirrec.fit = Mock()
        dirrec.curved = Mock()
        dirrec.direct.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        dirrec.fit.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        dirrec.curved.reconstruct_common.return_value = (sentinel.theta, sentinel.phi)
        coincidence_0 = []
        coincidence_1 = [[sentinel.station_number, {'timestamp': 1}]]
        coincidence_2 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event]]
        coincidence_3 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event]]
        coincidence_4 = [[sentinel.station_number, {'timestamp': 1}], [1, sentinel.event],
                         [2, sentinel.event], [3, sentinel.event]]

        # To few detection points, no reconstruction
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_0)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

        # To few detection points, no reconstruction
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_1)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 1)

        # To few eligible detection points, no reconstruction
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_2, station_numbers=[1])
        dirrec.fit.reconstruct_common.assert_not_called()
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 1)

        # Two stations with four detection points, fit reconstruction
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_2)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.fit.reconstruct_common.assert_called_once_with(
            [0.] * 4, [sentinel.x] * 4, [sentinel.y] * 4, [sentinel.z] * 4, {})
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 1)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 2)

        # Four events with eight detection points and initial core,
        # curved reconstruction
        initial = {'core_x': sentinel.core_x, 'core_y': sentinel.core_y}
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4, initial=initial)
        cluster.set_timestamp.assert_called_with(1)
        self.assertEqual(dirrec.curved.reconstruct_common.call_count, 1)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 4)

        mock_arrival_times.return_value = [0., nan, nan, nan]

        # Three stations with three detection points, direct reconstruction
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_3)
        cluster.set_timestamp.assert_called_with(1)
        dirrec.direct.reconstruct_common.assert_called_once_with(
            [0.] * 3, [sentinel.x] * 3, [sentinel.y] * 3, [sentinel.z] * 3, {})
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        self.assertEqual(theta, sentinel.theta)
        self.assertEqual(phi, sentinel.phi)
        self.assertEqual(len(nums), 3)

        mock_arrival_times.return_value = [nan] * 4

        # Four events, no valid station arrival times
        theta, phi, nums = dirrec.reconstruct_coincidence(coincidence_4)
        cluster.set_timestamp.assert_called_with(1)
        self.assertEqual(dirrec.curved.reconstruct_common.call_count, 1)
        self.assertEqual(dirrec.fit.reconstruct_common.call_count, 1)
        self.assertEqual(dirrec.direct.reconstruct_common.call_count, 1)
        self.assertTrue(isnan(theta))
        self.assertTrue(isnan(phi))
        self.assertEqual(len(nums), 0)

    @patch.object(direction_reconstruction.CoincidenceDirectionReconstructionDetectors, 'reconstruct_coincidence')
    def test_reconstruct_coincidences(self, mock_reconstruct_coincidence):
        mock_reconstruct_coincidence.return_value = [sentinel.theta, sentinel.phi, sentinel.nums]
        self.assertEqual(self.dirrec.reconstruct_coincidences([sentinel.coincidence, sentinel.coincidence],
                                                              sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((sentinel.theta, sentinel.theta), (sentinel.phi, sentinel.phi), (sentinel.nums, sentinel.nums)))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)
        mock_reconstruct_coincidence.assert_called_with(sentinel.coincidence, sentinel.station_numbers, sentinel.offsets, None)
        self.assertEqual(self.dirrec.reconstruct_coincidences([], sentinel.station_numbers, sentinel.offsets, progress=False),
                         ((), (), ()))
        self.assertEqual(mock_reconstruct_coincidence.call_count, 2)


class BaseAlgorithm(object):

    """Use this class to check the different algorithms

    This provides a shortcut to call the reconstruct_common method.

    """

    def call_reconstruct(self, t, x, y, z, initial=None):
        return self.algorithm.reconstruct_common(t, x, y, z, initial)


class FlatAlgorithm(BaseAlgorithm):

    """Use this class to test algorithms for flat shower fronts.

    They should give similar results and errors in some cases.
    These tests use three detections at same height.

    """

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
        self.assertTrue(-pi <= phi < pi)

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


class DirectAlgorithm(FlatAlgorithm):

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
            self.assertTrue(-pi <= phi < pi)


class AltitudeAlgorithm(FlatAlgorithm):

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


class MultiAlgorithm(FlatAlgorithm):

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


class CurvedAlgorithm(BaseAlgorithm):

    """Check some algorithms supporting a curved shower front.

    They should give similar results and errors in some cases.

    """

    def test_curved_shower(self):
        """Simple curved shower on three detectors."""

        t = (0., 0., 10.)
        x = (0., 100., 50.)
        y = (0., 0., 100.)
        z = (0., 0., 0.)
        init = {'core_x': 50, 'core_y': 0}

        theta, phi = self.call_reconstruct(t, x, y, z, initial=init)

        self.assertAlmostEqual(theta, 0., 4)
        self.assertTrue(-pi <= phi < pi)


class CurvedAltitudeAlgorithm(CurvedAlgorithm):

    """Check algorithms for curved fronts and stations at different altitudes.

    They should give similar results and errors in some cases.

    """

    def test_curved_shower_on_stations_with_altitude(self):
        """Simple curved shower on three stations at different altitudes."""

        c = 0.299792458

        z = (10, 0., 40.)
        t = (-z[0] / c, 0., 10. - z[2] / c)
        x = (0., 100., 50.)
        y = (0., 0., 100.)
        init = {'core_x': 50, 'core_y': 0}

        theta, phi = self.call_reconstruct(t, x, y, z, initial=init)

        self.assertAlmostEqual(theta, 0., 4)
        self.assertTrue(-pi <= phi < pi)


class DirectAlgorithmTest(unittest.TestCase, DirectAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithm()


class DirectAlgorithmCartesianTest(unittest.TestCase, DirectAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian()


class DirectAlgorithmCartesian3DTest(unittest.TestCase,
                                     DirectAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.DirectAlgorithmCartesian3D()


class FitAlgorithm3DTest(unittest.TestCase, MultiAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.FitAlgorithm3D()


class RegressionAlgorithmTest(unittest.TestCase, MultiAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.RegressionAlgorithm()


class RegressionAlgorithm3DTest(unittest.TestCase, MultiAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.RegressionAlgorithm3D()


class CurvedRegressionAlgorithmTest(unittest.TestCase, CurvedAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.CurvedRegressionAlgorithm()
        self.algorithm.front = ConeFront()


class CurvedRegressionAlgorithm3DTest(unittest.TestCase,
                                      CurvedAltitudeAlgorithm):

    def setUp(self):
        self.algorithm = direction_reconstruction.CurvedRegressionAlgorithm3D()
        self.algorithm.front = ConeFront()


if __name__ == '__main__':
    unittest.main()
