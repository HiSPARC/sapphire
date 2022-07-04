import os
import unittest

from unittest.mock import MagicMock, patch, sentinel

import tables

from sapphire.analysis import reconstructions

TEST_DATA_FILE = '../simulations/test_data/groundparticles_sim.h5'


class ReconstructESDEventsTest(unittest.TestCase):

    def setUp(self):
        self.data = MagicMock()
        self.station = MagicMock(spec=reconstructions.Station)
        self.rec = reconstructions.ReconstructESDEvents(
            self.data, sentinel.station_group, self.station,
            overwrite=sentinel.overwrite, progress=sentinel.progress,
            verbose=False, destination=sentinel.destination)

    def test_init(self):
        rec = self.rec
        self.assertEqual(rec.data, self.data)

        self.assertEqual(rec.station_group, self.data.get_node.return_value)
        self.data.get_node.assert_called_once_with(sentinel.station_group)
        self.assertEqual(rec.events, self.data.get_node.return_value.events)

        self.assertEqual(rec.overwrite, sentinel.overwrite)
        self.assertEqual(rec.progress, sentinel.progress)
        self.assertFalse(rec.verbose)
        self.assertEqual(rec.destination, sentinel.destination)
        self.assertEqual(rec.offsets, [0.] * 4)

        self.assertEqual(rec.station, self.station)

        self.assertTrue(isinstance(rec.direction, reconstructions.EventDirectionReconstruction))
        self.assertTrue(isinstance(rec.core, reconstructions.EventCoreReconstruction))

        self.assertEqual(rec.theta, [])
        self.assertEqual(rec.phi, [])
        self.assertEqual(rec.detector_ids, [])
        self.assertEqual(rec.core_x, [])
        self.assertEqual(rec.core_y, [])

    def test_reconstruct_directions(self):
        self.rec.direction = MagicMock()
        self.rec.direction.reconstruct_events.return_value = (sentinel.theta, sentinel.phi, sentinel.ids)
        self.rec.reconstruct_directions()
        self.rec.direction.reconstruct_events.assert_called_once_with(
            self.rec.events, None, self.rec.offsets, self.rec.progress, [])
        self.assertEqual(self.rec.theta, sentinel.theta)
        self.assertEqual(self.rec.phi, sentinel.phi)
        self.assertEqual(self.rec.detector_ids, sentinel.ids)

        self.rec.reconstruct_directions(sentinel.detector_ids)
        self.rec.direction.reconstruct_events.assert_called_with(
            self.rec.events, sentinel.detector_ids, self.rec.offsets, self.rec.progress, [])

    def test_reconstruct_cores(self):
        self.rec.core = MagicMock()
        self.rec.core.reconstruct_events.return_value = (sentinel.core_x, sentinel.core_y)
        self.rec.reconstruct_cores()
        self.rec.core.reconstruct_events.assert_called_once_with(
            self.rec.events, None, self.rec.progress, [])
        self.assertEqual(self.rec.core_x, sentinel.core_x)
        self.assertEqual(self.rec.core_y, sentinel.core_y)

        self.rec.reconstruct_cores(sentinel.detector_ids)
        self.rec.core.reconstruct_events.assert_called_with(
            self.rec.events, sentinel.detector_ids, self.rec.progress, [])

    def test_prepare_output(self):
        self.rec.events = MagicMock()
        self.rec.events.nrows = sentinel.nrows
        self.rec.prepare_output()
        self.data.create_table.assert_called_once_with(
            self.rec.station_group, sentinel.destination,
            reconstructions.ReconstructedEvent, expectedrows=sentinel.nrows)
        self.assertEqual(self.rec.reconstructions, self.data.create_table.return_value)
        self.assertEqual(self.rec.reconstructions._v_attrs.station, self.station)

    def test_prepare_output_existing(self):
        self.rec.events = MagicMock()
        self.rec.events.nrows = sentinel.nrows
        self.rec.station_group = [sentinel.destination]

        # Overwrite existing
        self.rec.overwrite = True
        self.rec.prepare_output()
        self.data.remove_node.assert_called_once_with(
            self.rec.station_group, sentinel.destination, recursive=True)
        self.data.create_table.assert_called_with(
            self.rec.station_group, sentinel.destination,
            reconstructions.ReconstructedEvent, expectedrows=sentinel.nrows)
        self.assertEqual(self.rec.reconstructions, self.data.create_table.return_value)
        self.assertEqual(self.rec.reconstructions._v_attrs.station, self.station)

        # Raise exception if table already exists
        self.rec.overwrite = False
        self.assertRaises(RuntimeError, self.rec.prepare_output)

    @patch.object(reconstructions.api, 'Station')
    @patch.object(reconstructions, 'determine_detector_timing_offsets')
    def test_get_detector_offsets(self, mock_determine_detctor_timing_offets, mock_station):
        mock_station.return_value = sentinel.station
        self.rec.events = sentinel.events
        self.rec.station.detectors = [None, None]

        # no offsets in station object no station_number ->
        #  determine offsets from events
        self.rec.get_detector_offsets()
        mock_determine_detctor_timing_offets.assert_called_with(
            sentinel.events, self.station)

        # no offsets in station object and station number -> api.Station
        self.rec.station_number = sentinel.station
        self.rec.get_detector_offsets()
        self.assertEqual(self.rec.offsets, sentinel.station)

        # offsets from cluster object (stored by simulation)
        detector = MagicMock(offset=sentinel.offset)
        self.rec.station = MagicMock(number=sentinel.number, detectors=[detector, detector])
        self.rec.get_detector_offsets()
        self.assertEqual(self.rec.offsets, [sentinel.offset, sentinel.offset])

    def test__store_reconstruction(self):
        event = MagicMock()
        # _store_reconstruction calls  min(event['n1'], ...).
        # but MagicMock is unordered in python 3!
        # Mock a dict that always returns 42.
        event.__getitem__.side_effect = lambda x: 42.
        self.rec.reconstructions = MagicMock()
        self.rec._store_reconstruction(event, sentinel.core_x, sentinel.core_y,
                                       sentinel.theta, sentinel.phi, [1, 3, 4])
        self.rec.reconstructions.row.append.assert_called_once_with()


class ReconstructESDEventsFromSourceTest(ReconstructESDEventsTest):

    def setUp(self):
        self.data = MagicMock()
        self.dest_data = MagicMock()
        self.station = MagicMock(spec=reconstructions.Station)
        self.rec = reconstructions.ReconstructESDEventsFromSource(
            self.data, self.dest_data, sentinel.station_group,
            sentinel.dest_group, self.station, overwrite=sentinel.overwrite,
            progress=sentinel.progress, verbose=False,
            destination=sentinel.destination)

    @unittest.skip('WIP')
    def test_prepare_output(self):
        pass

    @unittest.skip('WIP')
    def test_prepare_output_existing(self):
        pass


class ReconstructSimulatedEventsTest(unittest.TestCase):

    def setUp(self):
        self.data = None

    def tearDown(self):
        if isinstance(self.data, tables.file.File):
            self.data.close()

    def test_station_is_object(self):
        self.data = MagicMock()
        station = MagicMock(spec=reconstructions.Station)
        rec = reconstructions.ReconstructSimulatedEvents(
            self.data, sentinel.station_group, station,
            overwrite=sentinel.overwrite, progress=sentinel.progress,
            verbose=False, destination=sentinel.destination)
        self.assertEqual(rec.station, station)

    def test_read_object_from_hdf5(self):
        fn = self.get_testdata_path(TEST_DATA_FILE)
        self.data = tables.open_file(fn, 'r')
        station_group = '/cluster_simulations/station_0'
        rec = reconstructions.ReconstructSimulatedEvents(
            self.data, station_group, 0)

        # isinstance does not work on classes that are read from pickles.
        self.assertEqual(rec.station.station_id, 0)

        with self.assertRaises(RuntimeError):
            rec = reconstructions.ReconstructSimulatedEvents(
                self.data, station_group, -999)

    def get_testdata_path(self, fn):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, fn)


class ReconstructESDCoincidencesTest(unittest.TestCase):

    @patch.object(reconstructions, 'CoincidenceQuery')
    def setUp(self, mock_cq):
        self.data = MagicMock()
        self.cluster = MagicMock()
        self.cq = mock_cq
        self.rec = reconstructions.ReconstructESDCoincidences(
            self.data, sentinel.coin_group, overwrite=sentinel.overwrite,
            progress=sentinel.progress, verbose=False,
            destination=sentinel.destination, cluster=self.cluster)

    def test_init(self):
        rec = self.rec
        self.assertEqual(rec.data, self.data)

        self.assertEqual(rec.coincidences_group, self.data.get_node.return_value)
        self.data.get_node.assert_called_once_with(sentinel.coin_group)
        self.assertEqual(rec.coincidences, self.data.get_node.return_value.coincidences)

        self.assertEqual(rec.overwrite, sentinel.overwrite)
        self.assertEqual(rec.progress, sentinel.progress)
        self.assertFalse(rec.verbose)
        self.assertEqual(rec.destination, sentinel.destination)
        self.assertEqual(rec.offsets, {})

        self.cq.assert_called_once_with(self.data, rec.coincidences_group)
        self.assertEqual(self.rec.cq, self.cq.return_value)

        self.assertEqual(rec.cluster, self.cluster)

        self.assertTrue(isinstance(rec.direction, reconstructions.CoincidenceDirectionReconstruction))
        self.assertTrue(isinstance(rec.core, reconstructions.CoincidenceCoreReconstruction))

        self.assertEqual(rec.theta, [])
        self.assertEqual(rec.phi, [])
        self.assertEqual(rec.station_numbers, [])
        self.assertEqual(rec.core_x, [])
        self.assertEqual(rec.core_y, [])

    @patch.object(reconstructions.api, 'Station')
    def test_get_station_timing_offsets(self, mock_station):
        mock_station.return_value = sentinel.station
        station = MagicMock(number=sentinel.number, spec=['number'])
        self.rec.cluster.stations = [station]
        self.rec.get_station_timing_offsets()
        self.assertEqual(list(self.rec.offsets.keys()), [sentinel.number])
        self.assertEqual(list(self.rec.offsets.values()), [sentinel.station])

        detector = MagicMock(offset=1.)
        station = MagicMock(number=sentinel.number, gps_offset=2., detectors=[detector],
                            spec=['gps_offset', 'number'])
        self.rec.cluster.stations = [station]
        self.rec.get_station_timing_offsets()
        self.assertEqual(list(self.rec.offsets.keys()), [sentinel.number])
        self.assertEqual(list(self.rec.offsets.values()), [[3.]])

    def test_reconstruct_directions(self):
        self.rec.coincidences = MagicMock()
        self.rec.coincidences.nrows = 1
        self.rec.direction = MagicMock()
        self.rec.direction.reconstruct_coincidences.return_value = (sentinel.theta, sentinel.phi, sentinel.nums)
        self.rec.reconstruct_directions()
        self.rec.direction.reconstruct_coincidences.assert_called_once_with(
            self.rec.cq.all_events.return_value, None, self.rec.offsets, progress=False, initials=[])
        self.assertEqual(self.rec.theta, sentinel.theta)
        self.assertEqual(self.rec.phi, sentinel.phi)
        self.assertEqual(self.rec.station_numbers, sentinel.nums)

        self.rec.reconstruct_directions(sentinel.nums)
        self.rec.direction.reconstruct_coincidences.assert_called_with(
            self.rec.cq.all_events.return_value, sentinel.nums, self.rec.offsets, progress=False, initials=[])

    def test_reconstruct_cores(self):
        self.rec.coincidences = MagicMock()
        self.rec.coincidences.nrows = 1
        self.rec.core = MagicMock()
        self.rec.core.reconstruct_coincidences.return_value = (sentinel.core_x, sentinel.core_y)
        self.rec.reconstruct_cores()
        self.rec.core.reconstruct_coincidences.assert_called_once_with(
            self.rec.cq.all_events.return_value, None, progress=False, initials=[])
        self.assertEqual(self.rec.core_x, sentinel.core_x)
        self.assertEqual(self.rec.core_y, sentinel.core_y)

        self.rec.reconstruct_cores(sentinel.nums)
        self.rec.core.reconstruct_coincidences.assert_called_with(
            self.rec.cq.all_events.return_value, sentinel.nums, progress=False, initials=[])

    def test_prepare_output(self):
        self.rec.coincidences = MagicMock()
        self.rec.coincidences.nrows = sentinel.nrows
        self.cluster.stations.return_value = []
        self.rec.prepare_output()
        self.data.create_table.assert_called_once_with(
            self.rec.coincidences_group, sentinel.destination,
            reconstructions.ReconstructedCoincidence, expectedrows=sentinel.nrows)
        self.assertEqual(self.rec.reconstructions, self.data.create_table.return_value)
        self.assertEqual(self.rec.reconstructions._v_attrs.cluster, self.cluster)

    def test_prepare_output_existing(self):
        self.rec.coincidences = MagicMock()
        self.rec.coincidences.nrows = sentinel.nrows
        self.cluster.stations.return_value = []
        self.rec.coincidences_group = [sentinel.destination]

        # Overwrite existing
        self.rec.overwrite = True
        self.rec.prepare_output()
        self.data.remove_node.assert_called_once_with(
            self.rec.coincidences_group, sentinel.destination, recursive=True)
        self.data.create_table.assert_called_with(
            self.rec.coincidences_group, sentinel.destination,
            reconstructions.ReconstructedCoincidence, expectedrows=sentinel.nrows)
        self.assertEqual(self.rec.reconstructions, self.data.create_table.return_value)
        self.assertEqual(self.rec.reconstructions._v_attrs.cluster, self.cluster)

        # Raise exception if table already exists
        self.rec.overwrite = False
        self.assertRaises(RuntimeError, self.rec.prepare_output)

    @patch.object(reconstructions, 'ReconstructedCoincidence')
    def test_prepare_output_columns(self, mock_description):
        self.rec.coincidences = MagicMock()
        self.rec.coincidences.nrows = sentinel.nrows
        station = MagicMock()
        station.number = 1
        self.rec.cluster.stations = [station]

        self.rec.prepare_output()
        mock_description.columns.update.assert_called_once_with(
            {'s1': reconstructions.tables.BoolCol(pos=26)})
        self.data.create_table.assert_called_with(
            self.rec.coincidences_group, sentinel.destination, mock_description,
            expectedrows=sentinel.nrows)

    def test__store_reconstruction(self):
        coin = MagicMock()
        self.rec.reconstructions = MagicMock()
        self.rec._store_reconstruction(coin, sentinel.core_x, sentinel.core_y,
                                       sentinel.theta, sentinel.phi, [2, 3, 4])
        self.rec.reconstructions.row.append.assert_called_once_with()


class ReconstructESDCoincidencesFromSourceTest(ReconstructESDCoincidencesTest):

    @patch.object(reconstructions, 'CoincidenceQuery')
    def setUp(self, mock_cq):
        self.data = MagicMock()
        self.dest_data = MagicMock()
        self.cluster = MagicMock()
        self.cq = mock_cq
        self.rec = reconstructions.ReconstructESDCoincidencesFromSource(
            self.data, self.dest_data, sentinel.coin_group,
            sentinel.dest_group, overwrite=sentinel.overwrite,
            progress=sentinel.progress, verbose=False,
            destination=sentinel.destination, cluster=self.cluster)

    @unittest.skip('WIP')
    def test_prepare_output(self):
        pass

    @unittest.skip('WIP')
    def test_prepare_output_existing(self):
        pass

    @unittest.skip('WIP')
    def test_prepare_output_columns(self):
        pass


class ReconstructSimulatedCoincidencesTest(unittest.TestCase):

    def setUp(self):
        self.data = MagicMock()

    def tearDown(self):
        if isinstance(self.data, tables.file.File):
            self.data.close()

    @patch.object(reconstructions, 'CoincidenceQuery')
    def test_cluster_is_object(self, mock_cq):
        cluster = MagicMock()
        rec = reconstructions.ReconstructSimulatedCoincidences(
            self.data, sentinel.coin_group, overwrite=sentinel.overwrite,
            progress=sentinel.progress, verbose=False,
            destination=sentinel.destination, cluster=cluster)
        self.assertEqual(rec.cluster, cluster)

    def test_read_object_from_hdf5(self):
        fn = self.get_testdata_path(TEST_DATA_FILE)
        self.data = tables.open_file(fn, 'r')
        rec = reconstructions.ReconstructSimulatedCoincidences(self.data)

        # isinstance does not work on classes that are read from pickles.
        self.assertEqual(rec.cluster.stations[0].station_id, 0)

    def get_testdata_path(self, fn):
        dir_path = os.path.dirname(__file__)
        return os.path.join(dir_path, fn)


if __name__ == '__main__':
    unittest.main()
