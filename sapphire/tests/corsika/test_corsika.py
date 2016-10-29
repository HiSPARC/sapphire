import unittest
import os.path
from math import pi

from sapphire import corsika


data_file_dir = os.path.dirname(__file__)
DATA_FILE = os.path.join(data_file_dir, 'test_data/1_2/DAT000000')


class CorsikaFileTests(unittest.TestCase):
    def setUp(self):
        self.file = corsika.reader.CorsikaFile(DATA_FILE)

    def tearDown(self):
        self.file.finish()

    def test_validate_file(self):
        """Verify that the data file is valid"""

        self.assertTrue(self.file.check())

    def test_run_header(self):
        """Verify that the Run header is properly read"""

        header = self.file.get_header()
        self.assertIsInstance(header, corsika.blocks.RunHeader)
        self.assertEqual(header.id, b'RUNH')
        self.assertAlmostEqual(header.version, 7.4, 4)
        for h in [10., 5000., 30000., 50000., 110000.]:
            t = header.height_to_thickness(h)
            self.assertAlmostEqual(header.thickness_to_height(t), h, 8)

    def test_run_end(self):
        """Verify that the Run end is properly read"""

        end = self.file.get_end()
        self.assertIsInstance(end, corsika.blocks.RunEnd)
        self.assertEqual(end.id, b'RUNE')
        self.assertEqual(end.n_events_processed, 1)

    def test_events(self):
        """Verify that the Events are properly read"""

        events = self.file.get_events()
        event = next(events)
        self.assertIsInstance(event, corsika.reader.CorsikaEvent)
        self.assertEqual(event.last_particle_index, 1086892)

    def test_event_header(self):
        """Verify that the Event header is properly read"""

        events = self.file.get_events()
        event = next(events)
        header = event.get_header()
        self.assertIsInstance(header, corsika.blocks.EventHeader)
        self.assertEqual(header.id, b'EVTH')
        self.assertEqual(corsika.particles.name(header.particle_id), 'proton')
        self.assertEqual(header.energy, 1e14)
        self.assertEqual(header.azimuth, -pi / 2.)
        self.assertEqual(header.zenith, 0.0)
        self.assertEqual(header.hadron_model_high, 'QGSJET')

    def test_event_end(self):
        """Verify that the Event end is properly read"""

        events = self.file.get_events()
        event = next(events)
        end = event.get_end()
        self.assertIsInstance(end, corsika.blocks.EventEnd)
        self.assertEqual(end.id, b'EVTE')
        self.assertEqual(end.n_muons_output, 1729)

    def test_particles(self):
        """Verify that the Particles are properly read"""

        events = self.file.get_events()
        event = next(events)
        particles = event.get_particles()
        particle = next(particles)
        self.assertIsInstance(particle, tuple)
        self.assertEqual(len(particle), 11)
        self.assertEqual(corsika.particles.name(int(particle[6])), 'muon_p')
        self.assertAlmostEqual(particle[3], -56.2846679688)
        self.assertAlmostEqual(particle[4], -172.535859375)
        self.assertAlmostEqual(particle[7], 181.484397728)
        particle = next(particles)
        self.assertEqual(corsika.particles.name(int(particle[6])), 'muon_m')


if __name__ == '__main__':
    unittest.main()
