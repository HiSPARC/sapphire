import unittest

from math import atan2, sqrt

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

    def test_validate_block_format(self):
        """Verify that the block format is logical"""

        self.assertEqual(
            (self.format.block_size - 2 * self.format.block_padding_size) / self.format.subblock_size,
            self.format.subblocks_per_block,
            msg=(
                f'The block format ({self.format.block_format}) and sub-block format '
                f'({self.format.subblock_format}) do not agree! block size is {self.format.block_size} '
                f'and sub-block size is {self.format.subblock_size}. Block size should'
                f' be {self.format.subblocks_per_block} times the sub-block size plus '
                'padding (usually 8 bytes).'
            ),
        )

    def test_validate_subblock_format(self):
        """Verify that the subblock format is logical"""

        self.assertEqual(
            self.format.subblock_size / self.format.particle_size,
            self.format.particles_per_subblock,
            msg=(
                f'The sub_block format ({self.format.subblock_format}) and particle format '
                f'({self.format.particle_format}) do not agree! sub-block size is '
                f'{self.format.subblock_size} and particle record size is '
                f'{self.format.particle_size}. Sub-block size should be '
                f'{self.format.particles_per_subblock} times the particle record size.'
            ),
        )

    def test_validate_particle_format(self):
        """Verify that the particle format is correct"""

        self.assertEqual(
            self.format.particle_format,
            '7f',
            msg=(f'The particle format ({self.format.particle_format}) is incorrect.'),
        )


class CorsikaBlocksThinTests(CorsikaBlocksTests):
    def setUp(self):
        self.format = blocks.FormatThin()

    def test_validate_particle_format(self):
        """Verify that the particle format is correct"""

        self.assertEqual(
            self.format.particle_format,
            '8f',
            msg=(f'The thinned particle format ({self.format.particle_format}) is incorrect.'),
        )


class ParticleDataTests(unittest.TestCase):
    def setUp(self):
        # Input
        particle_id = 1000
        p_x = 2.0  # GeV
        p_y = 1.0  # GeV
        p_z = 10.0  # GeV
        x = 300.0  # cm
        y = 400.0  # cm
        t = 12345678.0  # ns

        self.subblock = (particle_id, p_x, p_y, p_z, x, y, t)

        # Output
        p_x *= 1e9  # eV
        p_y *= 1e9
        p_z *= 1e9
        x *= 1e-2  # m
        y *= 1e-2
        r = sqrt(x**2 + y**2)
        phi = atan2(x, -y)

        self.result = (p_x, p_y, -p_z, -y, x, t, particle_id / 1000, r, particle_id / 10 % 100, particle_id % 10, phi)

    def test_particle_data(self):
        """Verify conversion of particle information by particle_data()"""

        self.assertAlmostEqual(blocks.particle_data(self.subblock), self.result)

    @unittest.skipUnless(numba_available, 'Numba required')
    def test_numba_jit(self):
        """Verify particle_data() with numba JIT disabled"""

        self.assertTrue(hasattr(blocks.particle_data, '__numba__'))
        old_value = numba.config.DISABLE_JIT
        numba.config.DISABLE_JIT = 1
        self.assertAlmostEqual(blocks.particle_data(self.subblock), self.result)
        numba.config.DISABLE_JIT = old_value


class ParticleDataThinTests(ParticleDataTests):
    def setUp(self):
        super().setUp()

        # Input
        weight = 9.0
        self.subblock = self.subblock + (weight,)

        # Output
        self.result = self.result + (weight,)

    def test_particle_data(self):
        """Verify conversion of particle information by particle_data()"""

        self.assertAlmostEqual(blocks.particle_data_thin(self.subblock), self.result)

    @unittest.skipUnless(numba_available, 'Numba required')
    def test_numba_jit(self):
        """Verify particle_data() with numba JIT disabled"""

        self.assertTrue(hasattr(blocks.particle_data_thin, '__numba__'))
        old_value = numba.config.DISABLE_JIT
        numba.config.DISABLE_JIT = 1
        self.assertAlmostEqual(blocks.particle_data_thin(self.subblock), self.result)
        numba.config.DISABLE_JIT = old_value
