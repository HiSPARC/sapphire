import unittest

from mock import sentinel, MagicMock

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
        self.rec.prepare_output()
        self.data.create_table.assert_called_once_with(
            self.rec.station_group, sentinel.destination, reconstructions.ReconstructedEvent)
        self.assertEqual(self.rec.reconstructions, self.data.create_table.return_value)
        self.assertEqual(self.rec.reconstructions._v_attrs.station, self.station)

    def test_prepare_output_existing(self):
        self.rec.station_group = [sentinel.destination]

        # Overwrite existing
        self.rec.overwrite = True
        self.rec.prepare_output()
        self.data.remove_node.assert_called_once_with(
            self.rec.station_group, sentinel.destination, recursive=True)
        self.data.create_table.assert_called_with(
            self.rec.station_group, sentinel.destination, reconstructions.ReconstructedEvent)
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


if __name__ == '__main__':
    unittest.main()
