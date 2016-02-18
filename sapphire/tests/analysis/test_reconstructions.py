import unittest

from mock import sentinel, MagicMock, patch

from sapphire.analysis import reconstructions


class ReconstructESDEventsTest(unittest.TestCase):

    def setUp(self):
        self.data = MagicMock()
        self.station = MagicMock(spec=reconstructions.Station)
        self.rec = reconstructions.ReconstructESDEvents(
            self.data, sentinel.station_group, self.station,
            overwrite=sentinel.overwrite, progress=sentinel.progress,
            destination=sentinel.destination)

    def test_init(self):
        rec = self.rec
        self.assertEqual(rec.data, self.data)

        self.assertEqual(rec.station_group, self.data.get_node.return_value)
        self.data.get_node.assert_called_once_with(sentinel.station_group)
        self.assertEqual(rec.events, self.data.get_node.return_value.events)

        self.assertEqual(rec.overwrite, sentinel.overwrite)
        self.assertEqual(rec.progress, sentinel.progress)
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
            self.rec.events, None, self.rec.offsets, self.rec.progress)
        self.assertEqual(self.rec.theta, sentinel.theta)
        self.assertEqual(self.rec.phi, sentinel.phi)
        self.assertEqual(self.rec.detector_ids, sentinel.ids)

        self.rec.reconstruct_directions(sentinel.detector_ids)
        self.rec.direction.reconstruct_events.assert_called_with(
            self.rec.events, sentinel.detector_ids, self.rec.offsets, self.rec.progress)

    def test_reconstruct_cores(self):
        self.rec.core = MagicMock()
        self.rec.core.reconstruct_events.return_value = (sentinel.core_x, sentinel.core_y)
        self.rec.reconstruct_cores()
        self.rec.core.reconstruct_events.assert_called_once_with(self.rec.events, None, self.rec.progress)
        self.assertEqual(self.rec.core_x, sentinel.core_x)
        self.assertEqual(self.rec.core_y, sentinel.core_y)

        self.rec.reconstruct_cores(sentinel.detector_ids)
        self.rec.core.reconstruct_events.assert_called_with(
            self.rec.events, sentinel.detector_ids, self.rec.progress)

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

    def test__store_reconstruction(self):
        event = MagicMock()
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
            progress=sentinel.progress, destination=sentinel.destination)

    @unittest.skip('WIP')
    def test_prepare_output(self):
        pass

    @unittest.skip('WIP')
    def test_prepare_output_existing(self):
        pass


class ReconstructESDCoincidencesTest(unittest.TestCase):

    @patch.object(reconstructions, 'CoincidenceQuery')
    def setUp(self, mock_cq):
        self.data = MagicMock()
        self.cluster = MagicMock()
        self.cq = mock_cq
        self.rec = reconstructions.ReconstructESDCoincidences(
            self.data, sentinel.coin_group, overwrite=sentinel.overwrite,
            progress=sentinel.progress, destination=sentinel.destination,
            cluster=self.cluster)

    def test_init(self):
        rec = self.rec
        self.assertEqual(rec.data, self.data)

        self.assertEqual(rec.coincidences_group, self.data.get_node.return_value)
        self.data.get_node.assert_called_once_with(sentinel.coin_group)
        self.assertEqual(rec.coincidences, self.data.get_node.return_value.coincidences)

        self.assertEqual(rec.overwrite, sentinel.overwrite)
        self.assertEqual(rec.progress, sentinel.progress)
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

    def test_reconstruct_directions(self):
        self.rec.coincidences = MagicMock()
        self.rec.coincidences.nrows = 1
        self.rec.direction = MagicMock()
        self.rec.direction.reconstruct_coincidences.return_value = (sentinel.theta, sentinel.phi, sentinel.nums)
        self.rec.reconstruct_directions()
        self.rec.direction.reconstruct_coincidences.assert_called_once_with(
            self.rec.cq.all_events.return_value, None, self.rec.offsets, progress=False)
        self.assertEqual(self.rec.theta, sentinel.theta)
        self.assertEqual(self.rec.phi, sentinel.phi)
        self.assertEqual(self.rec.station_numbers, sentinel.nums)

        self.rec.reconstruct_directions(sentinel.nums)
        self.rec.direction.reconstruct_coincidences.assert_called_with(
            self.rec.cq.all_events.return_value, sentinel.nums, self.rec.offsets, progress=False)

    def test_reconstruct_cores(self):
        self.rec.coincidences = MagicMock()
        self.rec.coincidences.nrows = 1
        self.rec.core = MagicMock()
        self.rec.core.reconstruct_coincidences.return_value = (sentinel.core_x, sentinel.core_y)
        self.rec.reconstruct_cores()
        self.rec.core.reconstruct_coincidences.assert_called_once_with(
            self.rec.cq.all_events.return_value, None, progress=False)
        self.assertEqual(self.rec.core_x, sentinel.core_x)
        self.assertEqual(self.rec.core_y, sentinel.core_y)

        self.rec.reconstruct_cores(sentinel.nums)
        self.rec.core.reconstruct_coincidences.assert_called_with(
            self.rec.cq.all_events.return_value, sentinel.nums, progress=False)

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
            progress=sentinel.progress, destination=sentinel.destination,
            cluster=self.cluster)

    @unittest.skip('WIP')
    def test_prepare_output(self):
        pass

    @unittest.skip('WIP')
    def test_prepare_output_existing(self):
        pass

    @unittest.skip('WIP')
    def test_prepare_output_columns(self):
        pass


if __name__ == '__main__':
    unittest.main()
