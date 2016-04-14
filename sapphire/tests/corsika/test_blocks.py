import unittest
from mock import patch
from math import sqrt, atan2

from sapphire.corsika import blocks

try:
    import numba
    numba.__version__  # stop flake8 from complaining about unused module
except ImportError:
    numba_available = False
else:
    numba_available = True


class CorsikaBlocksTests(unittest.TestCase):
    def setUp(self):
        self.format = blocks.Format()

    def tearDown(self):
        pass

    def test_validate_block_format(self):
        """Verify that the block format is logical"""

        self.assertEqual((self.format.block_size - 2 * self.format.block_padding_size) / self.format.subblock_size,
                         self.format.subblocks_per_block,
                         msg=('The block format ({block}) and sub-block format '
                              '({sub_block}) do not agree! block size is {block_size} '
                              'and sub-block size is {sub_block_size}. Block size should'
                              ' be {subblocks_per_block} times the sub-block size plus '
                              'padding (usually 8 bytes).'
                              .format(block=self.format.block_format,
                                      sub_block=self.format.subblock_format,
                                      block_size=self.format.block_size,
                                      sub_block_size=self.format.subblock_size,
                                      subblocks_per_block=self.format.subblocks_per_block)))

    def test_validate_subblock_format(self):
        """Verify that the subblock format is logical"""

        self.assertEqual(self.format.subblock_size / self.format.particle_size,
                         self.format.particles_per_subblock,
                         msg=('The sub_block format ({sub_block}) and particle format '
                              '({particle}) do not agree! sub-block size is '
                              '{sub_block_size} and particle record size is '
                              '{particle_size}. Sub-block size should be '
                              '{particles_per_subblock} times the particle record size.'
                              .format(sub_block=self.format.subblock_format,
                                      particle=self.format.particle_format,
                                      sub_block_size=self.format.subblock_size,
                                      particle_size=self.format.particle_size,
                                      particles_per_subblock=self.format.particles_per_subblock)))

    def test_validate_particle_format(self):
        """Verify that the particle format is correct"""

        self.assertEqual(self.format.particle_format, '7f',
                         msg=('The particle format ({particle}) is incorrect.'
                              .format(particle=self.format.particle_format)))


class CorsikaBlocksThinTests(CorsikaBlocksTests):
    def setUp(self):
        self.format = blocks.FormatThin()

    def test_validate_particle_format(self):
        """Verify that the particle format is correct"""

        self.assertEqual(self.format.particle_format, '8f',
                         msg=('The thinned particle format ({particle}) is incorrect.'
                              .format(particle=self.format.particle_format)))


class ParticleDataTests(unittest.TestCase):
    def setUp(self):
        id = 1000
        p_x = 2.   # GeV
        p_y = 1.   # GeV
        p_z = 10.  # GeV
        x = 300.   # cm
        y = 400.   # cm
        t = 12345678.  # ns

        self.subblock = (id, p_x, p_y, p_z, x, y, t)

        p_x *= 1e9  # eV
        p_y *= 1e9
        p_z *= 1e9
        x *= 1e-2   # m
        y *= 1e-2

        self.result = (p_x, p_y, -p_z, -y, x, t, id / 1000,
                       sqrt(x ** 2 + y ** 2), id / 10 % 100,
                       id % 10, atan2(x, -y))

    def test_particle_data(self):
        """ verify conversion of particle information by particle_data() """
        self.assertAlmostEqual(blocks.particle_data(self.subblock), self.result)


class ParticleDataNumbaTests(ParticleDataTests):
    @unittest.skipUnless(numba_available, "Numba required")
    def test_numba_jit(self):
        """ verify particle_data() when compiled by numba """
        self.assertTrue(hasattr(blocks.particle_data, '__numba__'))
        self.assertAlmostEqual(blocks.particle_data(self.subblock), self.result)

    @unittest.skipUnless(numba_available, "Numba required")
    @patch('numba.jit', lambda x: x)
    def test_without_jit(self):
        """ test particla_data() without numba, while numba is installed """
        reload(blocks)
        self.assertFalse(hasattr(blocks.particle_data, '__numba__'))
        self.assertAlmostEqual(blocks.particle_data(self.subblock), self.result)

if __name__ == '__main__':
    unittest.main()
