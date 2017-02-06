import unittest

import numpy as np
from mock import patch

from sapphire.simulations import gammas


class GammasTest(unittest.TestCase):

    def test_compton_edge(self):
        # Compton edges of well known gammas sources
        # http://web.mit.edu/lululiu/Public/8.13/xray/TKA%20files/annihilation-Na.pdf
        # Na-22
        # 0.511 MeV: 0.344 MeV
        # 1.27 MeV: 1.06 MeV
        # http://www.spectrumtechniques.com/PDF/Compton%20Scattering%20Experiment%20by%20Prutchi.pdf
        # Cs-137
        # 0.662 MeV : 482 MeV
        # Co-60
        # 1.17 MeV: 0.96 MeV
        # 1.33 MeV: 1.12 MeV
        combinations = ((0.511, 0.340), (1.27, 1.06),
                        (0.662, 0.482),
                        (1.17, 0.96), (1.33, 1.12))

        for E, edge in combinations:
            self.assertAlmostEqual(gammas.compton_edge(E), edge, places=2)

    def test_compton_mean_free_path(self):
        # Relevant mean-free-paths in vinyltoluene scintillator
        # Values checked with: Jos Steijer, Nikhef internal note, 16 juni 2010, figure 3
        combinations = ((1., 32.), (10., 60.))

        for E, edge in combinations:
            self.assertAlmostEqual(gammas.compton_mean_free_path(E), edge, places=0)

    def test_pair_mean_free_path(self):
        # Relevant mean-free-paths in vinyltoluene scintillator
        # Values checked with: Jos Steijer, Nikhef internal note, 16 juni 2010, figure 5
        combinations = ((10, 249.), (1000., 62.))

        for E, edge in combinations:
            self.assertAlmostEqual(gammas.pair_mean_free_path(E), edge, places=0)

    @patch.object(np.random, 'random')
    def test_compton_energy_transfer(self, mock_random):
        # if random() return 1, energy should approach the kinematic maximum (compton edge)
        mock_random.return_value = 1.0
        for gamma_energy in [3., 10., 100.]:
            expected = gammas.compton_edge(gamma_energy)
            self.assertAlmostEqual(gammas.compton_energy_transfer(gamma_energy), expected, places=0)

        # if random() returns 0, energy should approach 0
        mock_random.return_value = 0.0
        for gamma_energy in [3., 10., 100.]:
            self.assertAlmostEqual(gammas.compton_energy_transfer(gamma_energy), 0.)

    def test_energy_transfer_cross_section(self):
        # The plot from github.com/tomkooij/lio-project/photons/check_sapphire_gammas.py
        # has been checked with Evans (1955) p. 693 figure 5.1
        # Peaks from this plot: (E [MeV], cross_section [barn])
        combinations = ((0.511, 1.6), (1.2, 0.52), (2.76, 0.22))

        barn = 1e-28  # m**2. Note that the figure in Evans is in centibarn!
        for E, cross_section in combinations:
            edge = gammas.compton_edge(E)
            self.assertAlmostEqual(gammas.energy_transfer_cross_section(E, edge) / barn, cross_section, places=1)

    def test_max_energy_transfer(self):
        self.assertAlmostEqual(gammas.max_energy_deposit_in_mips(0., 1.), gammas.MAX_E / gammas.MIP)
        self.assertAlmostEqual(gammas.max_energy_deposit_in_mips(0.5, 1.), 0.5 * gammas.MAX_E / gammas.MIP)
        self.assertAlmostEqual(gammas.max_energy_deposit_in_mips(1., 1.), 0.)

    @patch.object(gammas, 'compton_energy_transfer')
    @patch.object(gammas, 'pair_mean_free_path')
    @patch.object(gammas, 'compton_mean_free_path')
    def test_simulate_detector_mips_gammas_compton(self, mock_l_compton, mock_l_pair, mock_compton):
        # Force compton scattering
        mock_l_compton.return_value = 1e-3
        mock_l_pair.return_value = 1e50

        mock_compton.return_value = 1.
        p = np.array([10])
        theta = np.array([0.])

        mips = gammas.simulate_detector_mips_gammas(p, theta)
        mock_compton.assert_called_once_with(10. / 1e6)
        self.assertLessEqual(mips, gammas.MAX_E)

    @patch.object(gammas, 'compton_energy_transfer')
    @patch.object(gammas, 'pair_mean_free_path')
    @patch.object(gammas, 'compton_mean_free_path')
    def test_simulate_detector_mips_gammas_pair(self, mock_l_compton, mock_l_pair, mock_compton):
        # Force pair production
        mock_l_compton.return_value = 1e50
        mock_l_pair.return_value = 1e-3

        mock_compton.return_value = 42.
        E = np.array([10., 7.])  # MeV
        p = E * 1e6  # eV
        theta = np.array([0.])

        for _ in range(100):
            mips = gammas.simulate_detector_mips_gammas(p, theta)
            self.assertFalse(mock_compton.called)
            self.assertLessEqual(mips, gammas.MAX_E)

        # not enough energy for pair production
        E = np.array([0.5, 0.7])  # MeV
        p = E * 1e6  # eV
        theta = np.array([0., 0.])
        for _ in range(100):
            self.assertEqual(gammas.simulate_detector_mips_gammas(p, theta), 0)

    @patch('sapphire.simulations.gammas.expovariate')
    def test_simulate_detector_mips_no_interaction(self, mock_expovariate):
        p = np.array([10e6])
        theta = np.array([0.])

        # force no interaction
        mock_expovariate.side_effect = [1e6, 1e3]
        self.assertEqual(gammas.simulate_detector_mips_gammas(p, theta), 0)

        # no interaction because after projected depth
        mock_expovariate.side_effect = [4, 5]
        theta = np.array([1])
        self.assertEqual(gammas.simulate_detector_mips_gammas(p, theta), 0)

        # interactions are to late because of max depth
        mock_expovariate.side_effect = [120, 125]
        theta = np.array([1.555])  # projected depth would be 126 cm
        self.assertEqual(gammas.simulate_detector_mips_gammas(p, theta), 0)

        # no interaction with multiple inclined gammas each with to long
        # interaction depth.
        # test mostly to prevent accidental growing of depth in loop
        n = 30
        mock_expovariate.side_effect = [4, 5] * n
        p = np.array([10e6] * n)
        theta = np.array([1.] * n)  # projected depth would be 126 cm
        self.assertEqual(gammas.simulate_detector_mips_gammas(p, theta), 0)


if __name__ == '__main__':
    unittest.main()
