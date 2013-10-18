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

        self.assertTrue(self.file.Check())

    def test_run(self):
        """Verify that the Run is properly read"""

        header = self.file.get_header()
        self.assertIsInstance(header, corsika.blocks.RunHeader)
        self.assertEqual(header.fId, 'RUNH')
        self.assertAlmostEqual(header.fVersion, 7.4, 4)

        trailer = self.file.get_trailer()
        self.assertIsInstance(trailer, corsika.blocks.RunTrailer)
        self.assertEqual(trailer.fId, 'RUNE')

    def test_events(self):
        """Verify that the Events are properly read"""

        events = self.file.GetEvents()
        event = events.next()
        self.assertIsInstance(event, corsika.reader.CorsikaEvent)
        self.assertEqual(event.fLastParticle, 1086892)

    def test_event_header(self):
        """Verify that the Event header is properly read"""

        events = self.file.GetEvents()
        event = events.next()
        header = event.GetHeader()
        self.assertIsInstance(header, corsika.blocks.EventHeader)
        self.assertEqual(header.fEnergy, 1e14)

    def test_event_trailer(self):
        """Verify that the Event trailer is properly read"""

        events = self.file.GetEvents()
        event = events.next()
        trailer = event.GetTrailer()
        self.assertIsInstance(trailer, corsika.blocks.EventTrailer)
        self.assertEqual(trailer.fMuons, 1729)

    def test_particles(self):
        """Verify that the Particles are properly read"""

        events = self.file.GetEvents()
        event = events.next()
        particles = event.GetParticles()
        particle = particles.next()
        self.assertIsInstance(particle, corsika.blocks.ParticleData)
        self.assertEqual(particle.fParticle, corsika.particles.muon_p)
        particle = particles.next()
        self.assertEqual(particle.fParticle, corsika.particles.muon_m)


if __name__ == '__main__':
    unittest.main()
