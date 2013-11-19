import unittest

from sapphire.corsika import blocks


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


if __name__ == '__main__':
    unittest.main()
