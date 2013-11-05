import unittest

from sapphire import corsika


DATA_FILE = 'DAT000000'


class CorsikaFileTests(unittest.TestCase):
    def setUp(self):
        self.file = corsika.CorsikaFile(DATA_FILE)

    def tearDown(self):
        pass

    def test_validate_file(self):
        """Verify that the data file is valid"""

        self.assertTrue(self.file.check())

    def test_run(self):
        """Verify that the Run is properly read"""

        header = self.file.get_header()
        self.assertIsInstance(header, corsika.blocks.RunHeader)
        self.assertEqual(header.id, 'RUNH')
        self.assertAlmostEqual(header.version, 7.4, 4)

        end = self.file.get_end()
        self.assertIsInstance(end, corsika.blocks.RunEnd)
        self.assertEqual(end.id, 'RUNE')
        self.assertEqual(end.n_events_processed, 1)

    def test_events(self):
        """Verify that the Events are properly read"""

        events = self.file.get_events()
        event = events.next()
        self.assertIsInstance(event, corsika.reader.CorsikaEvent)
        self.assertEqual(event.last_particle_index, 1086892)

    def test_event_header(self):
        """Verify that the Event header is properly read"""

        events = self.file.get_events()
        event = events.next()
        header = event.get_header()
        self.assertIsInstance(header, corsika.blocks.EventHeader)
        self.assertEqual(header.id, 'EVTH')
        self.assertEqual(corsika.particles.id[header.particle_id], 'proton')
        self.assertEqual(header.energy, 1e14)
        self.assertEqual(header.hadron_model_high(), 'QGSJET')

    def test_event_end(self):
        """Verify that the Event end is properly read"""

        events = self.file.get_events()
        event = events.next()
        end = event.get_end()
        self.assertIsInstance(end, corsika.blocks.EventEnd)
        self.assertEqual(end.id, 'EVTE')
        self.assertEqual(end.n_muons_output, 1729)

    def test_particles(self):
        """Verify that the Particles are properly read"""

        events = self.file.get_events()
        event = events.next()
        particles = event.get_particles()
        particle = particles.next()
        self.assertIsInstance(particle, corsika.blocks.ParticleData)
        self.assertEqual(corsika.particles.id[particle.id], 'muon_p')
        self.assertAlmostEqual(particle.x, -172.535859375)
        self.assertAlmostEqual(particle.y, 56.2846679688)
        self.assertAlmostEqual(particle.r, 181.484397728)
        particle = particles.next()
        self.assertEqual(corsika.particles.id[particle.id], 'muon_m')


if __name__ == '__main__':
    unittest.main()
