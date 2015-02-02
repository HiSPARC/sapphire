import unittest

from sapphire.corsika import particles


class CorsikaParticlesTests(unittest.TestCase):

    def setUp(self):
        self.pid_name = [(1, 'gamma'),
                         (2, 'positron'),
                         (3, 'electron'),
                         (5, 'muon_p'),
                         (6, 'muon_m'),
                         (13, 'neutron'),
                         (14, 'proton'),
                         (201, 'deuteron'),
                         (301, 'tritium'),
                         (302, 'helium3'),
                         (402, 'alpha'),
                         (1206, 'carbon'),
                         (1407, 'nitrogen'),
                         (1608, 'oxygen'),
                         (2713, 'aluminium'),
                         (2814, 'silicon'),
                         (3216, 'sulfur'),
                         (5626, 'iron')]
        self.massless_atoms = [(909, 'fluorine'),
                               (3232, 'germanium'),
                               (9999, 'einsteinium')]
        self.atoms = [(1406, 'carbon14'),
                      (9999, 'einsteinium99')]

    def test_particle_IDs(self):
        """Verify that the correct names belong to each ID"""

        for id, name in self.pid_name:
            self.assertEqual(particles.ID[id], name)

    def test_conversion_functions(self):
        """Verify that the functions correctly convert back and forth"""

        for id, name in self.pid_name:
            self.assertEqual(particles.name(id), name)
            self.assertEqual(particles.particle_id(name), id)

        for id, name in self.massless_atoms:
            self.assertEqual(particles.particle_id(name), id)

        for id, name in self.atoms:
            self.assertEqual(particles.name(id), name)
            self.assertEqual(particles.particle_id(name), id)


if __name__ == '__main__':
    unittest.main()
