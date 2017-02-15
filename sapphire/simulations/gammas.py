"""Gammas

This module implements simulation of interactions of gamma photons with
a scintillator. Used in groundparticle simulations.

"""
from __future__ import division

from random import expovariate

import numpy as np

SCINTILLATOR_THICKNESS = 2.0  # cm
MAX_DEPTH = 112.  # longest straight path in scintillator in cm
ENERGY_LOSS = 2.0  # 2 MeV per cm
MAX_E = ENERGY_LOSS * SCINTILLATOR_THICKNESS
MIP = 3.38  # MeV

ELECTRON_REST_MASS_MeV = 0.5109989  # MeV


def compton_edge(gamma_energy):
    """Calculate Compton edge for a given photon energy

    W.R. Leo (1987) p.54

    :param gamma_energy: photon energy [MeV].
    :return: compton edge [MeV].

    """
    gamma = gamma_energy / ELECTRON_REST_MASS_MeV

    return gamma_energy * 2 * gamma / (1 + 2 * gamma)


def compton_energy_transfer(gamma_energy):
    """Calculate the energy transfer from photon to electron

    From the differential cross section the cumulative distribution
    is calculated. From this distribution a random energy transfer
    (within kinematic bounds) is returned.

    :param gamma_energy: photon energy [MeV].
    :return: transfered energy [MeV].

    """
    edge = compton_edge(gamma_energy)
    recoil_energies = np.linspace(0, edge, 1000)

    # electron energy distribution
    electron_energy = [energy_transfer_cross_section(gamma_energy,
                                                     recoil_energy)
                       for recoil_energy in recoil_energies]

    cumulative_energy = np.cumsum(electron_energy)

    normalised_energy_distribution = (cumulative_energy /
                                      cumulative_energy[-1])

    r = np.random.random()
    conversion_factor = normalised_energy_distribution.searchsorted(r) / 1000
    return compton_edge(gamma_energy) * conversion_factor


def energy_transfer_cross_section(gamma_energy, recoil_energy):
    """Differential cross section dsigma/dT

    Differential cross section for energy transfer from gamma
    to scattered electron in compton scattering.

    W.R. Leo (1987) p 54

    :param gamma_energy: photon energy [MeV].
    :param recoil_energy: electron recoil energy [MeV].

    """
    r_e = 2.82e-15  # classical electron radius [m]

    gamma = gamma_energy / ELECTRON_REST_MASS_MeV

    s = recoil_energy / gamma_energy

    return (np.pi * (r_e ** 2) / (ELECTRON_REST_MASS_MeV * gamma ** 2) *
            (2 + (s ** 2 / ((gamma ** 2) * ((1 - s) ** 2))) +
            (s / (1 - s)) * (s - 2 / gamma)))


def max_energy_deposit_in_mips(depth, scintillator_depth):
    """Maximum energy transfer from electron to scintillator

    Determine maximum energy transfer based on remaining scinitillator
    depth.

    Assumes scintillator depth is projected onto the direction
    of the incident particle (divided by cos(theta)).

    :param depth: depth at which the electron is produced [cm].
    :param scintillator_depth: total depth of the scintillator [cm].

    """
    return (scintillator_depth - depth) * MAX_E / (scintillator_depth * MIP)


def simulate_detector_mips_gammas(p, theta):
    """Simulate detection of gammas

    :param p: the momenta of the gammas as array, in eV.
    :param theta: angles of incidence of the gammas as array, in radians.
    :return: the simulated detector signal (in mips).

    """
    # p [eV] and E [MeV]
    energies = p / 1e6

    mips = 0
    for energy, angle in zip(energies, theta):
        # project depth onto direction of incident particle
        scintillator_depth = min(SCINTILLATOR_THICKNESS / np.cos(angle),
                                 MAX_DEPTH)

        # Calculate interaction point in units of scinitlator depth.
        # If depth > 1 there is no interaction.
        depth_compton = expovariate(1 / compton_mean_free_path(energy))
        depth_pair = expovariate(1 / pair_mean_free_path(energy))

        if ((depth_pair > scintillator_depth) &
                (depth_compton > scintillator_depth)):
            # no interaction
            continue

        # Interactions in scintillator
        elif depth_compton < depth_pair:
            # Compton scattering

            # kinetic energy transfered to electron by compton scattering
            energy_deposit = compton_energy_transfer(energy) / MIP
            max_deposit = max_energy_deposit_in_mips(depth_compton,
                                                     scintillator_depth)
            mips += min(max_deposit, energy_deposit)

        elif energy > 1.022:
            # Pair production: Two "electrons"

            # 1.022 MeV used for creation of two particles
            # all the rest is electron kinetic energy
            energy_deposit = (energy - 1.022) / MIP
            max_deposit = max_energy_deposit_in_mips(depth_pair,
                                                     scintillator_depth)
            mips += min(max_deposit, energy_deposit)

    return mips


def pair_mean_free_path(gamma_energy):
    """Mean free path pair production

    NIST XCOM database: http://www.nist.gov/pml/data/xcom/
    compound: C9H10
    pair production (total attenuation)

    table generated by @tomkooij/lio-project/photons/nist.py

    :param gamma_energy: photon energy [MeV].
    :return: mean free path [cm].

    """
    energy_path_pair_production = np.array([
        (4, 689.31), (5, 504.52), (6, 404.96),
        (7, 343.56), (8, 302.00), (9, 271.84),
        (10, 249.03), (11, 231.28), (12, 217.04),
        (13, 205.23), (14, 195.32), (15, 186.88),
        (16, 179.47), (18, 167.40), (20, 157.85),
        (22, 149.97), (24, 143.51), (26, 138.00),
        (28, 133.30), (30, 129.20), (40, 114.65),
        (50, 105.64), (60, 99.37), (80, 91.17),
        (100, 85.90), (150, 78.25), (200, 74.07),
        (300, 69.44), (400, 66.93), (500, 65.34),
        (600, 64.21), (800, 62.73), (1000, 61.82),
        (1500, 60.47), (2000, 59.72), (3000, 58.97),
        (4000, 58.53), (5000, 58.28), (6000, 58.09),
        (8000, 57.85), (10000, 57.70), (15000, 57.51),
        (20000, 57.41), (30000, 57.27), (40000, 57.21),
        (50000, 57.17), (60000, 57.13), (80000, 57.12),
        (100000, 57.08)])

    gamma_energies = energy_path_pair_production[:, 0]
    mean_free_paths = energy_path_pair_production[:, 1]

    idx = gamma_energies.searchsorted(gamma_energy, side='left')
    return mean_free_paths[idx]


def compton_mean_free_path(gamma_energy):
    """Mean free path compton scattering

    NIST XCOM database: http://www.nist.gov/pml/data/xcom/
    compound: C9H10
    compton scattering (incoherent scattering)

    table generated by @tomkooij/lio-project/photons/nist.py

    :param gamma_energy: photon energy [MeV].
    :return: mean free path [cm].

    """
    energy_path_compton_scattering = np.array([
        (4, 31.88), (5, 36.90), (6, 41.75),
        (7, 46.47), (8, 51.05), (9, 55.52),
        (10, 59.95), (11, 64.27), (12, 68.54),
        (13, 72.73), (14, 76.86), (15, 80.97),
        (16, 85.03), (18, 93.02), (20, 100.92),
        (22, 108.60), (24, 116.23), (26, 123.81),
        (28, 131.23), (30, 138.64), (40, 174.40),
        (50, 208.94), (60, 242.54), (80, 307.50),
        (100, 370.51), (150, 520.29), (200, 663.57),
        (300, 936.33), (400, 1195.46), (500, 1444.04),
        (600, 1686.34), (800, 2159.36), (1000, 2624.67),
        (1500, 3757.99), (2000, 4856.73), (3000, 6983.24),
        (4000, 9049.77), (5000, 11063.17), (6000, 13048.02),
        (8000, 16940.54), (10000, 20746.89), (15000, 30021.01),
        (20000, 39047.25), (30000, 56625.14), (40000, 73746.31),
        (50000, 90579.71), (60000, 107146.68), (80000, 139684.31),
        (100000, 171791.79)])

    gamma_energies = energy_path_compton_scattering[:, 0]
    mean_free_paths = energy_path_compton_scattering[:, 1]

    idx = gamma_energies.searchsorted(gamma_energy, side='left')
    return mean_free_paths[idx]
