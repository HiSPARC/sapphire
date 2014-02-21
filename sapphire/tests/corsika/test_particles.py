import math
import unittest

from sapphire.corsika import particles


class CorsikaParticlesTests(unittest.TestCase):

    def test_particle_ids(self):
        """Verify that the correct names belong to each id"""

        # Gamma
        self.assertEqual(particles.id[1], 'gamma')

        # Leptons
        self.assertEqual(particles.id[2], 'positron')
        self.assertEqual(particles.id[3], 'electron')
        self.assertEqual(particles.id[5], 'muon_p')
        self.assertEqual(particles.id[6], 'muon_m')

        # Hadrons
        self.assertEqual(particles.id[13], 'neutron')
        self.assertEqual(particles.id[14], 'proton')

        # Nucleons
        self.assertEqual(particles.id[201], 'deuteron')
        self.assertEqual(particles.id[301], 'tritium')
        self.assertEqual(particles.id[302], 'helium3')
        self.assertEqual(particles.id[402], 'alpha')
        self.assertEqual(particles.id[1206], 'carbon')
        self.assertEqual(particles.id[1407], 'nitrogen')
        self.assertEqual(particles.id[1608], 'oxygen')
        self.assertEqual(particles.id[2713], 'aluminium')
        self.assertEqual(particles.id[2814], 'silicon')
        self.assertEqual(particles.id[3216], 'sulfur')
        self.assertEqual(particles.id[5626], 'iron')


if __name__ == '__main__':
    unittest.main()
