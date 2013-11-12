"""
Classes corresponding to CORSIKA blocks and sub-blocks

The classes in this module correspond one-to-one with the sub-blocks
(and blocks) as specified in the CORSIKA users manual.

Author: Javier Gonzalez <jgonzalez@ik.fzk.de>
"""

import textwrap
import struct

import numpy

import units
import particles


# All sizes are in bytes

class Format(object):
    """
    Class containing the format information of the file
    as specified in Corsika user manual, Section 10.2.1.

    """
    def __init__(self):
        # in 32 bit, one field is a float
        self.field_size = struct.calcsize('f')

        # one block contains 5733 fields plus one header and one end field
        self.block_format = '5735f'
        self.block_size = struct.calcsize(self.block_format)
        self.block_padding_size = struct.calcsize('f')

        # Each block contains 21 sub-blocks
        # Each sub-block consists of 273 fields
        # the first of which _might_ be a string id.
        self.subblock_format = '4s272f'
        self.subblock_size = struct.calcsize(self.subblock_format)
        self.subblocks_per_block = 21

        # Each particle record sub-block contains a fixed
        # number of particle records
        # With the unthinned option, each of these is 7 fields long
        # for a total of 39 records per sub block
        self.particle_format = '7f'
        self.particle_size = struct.calcsize(self.particle_format)
        self.particles_per_subblock = 39


# From here on, things should not depend on the field size as everything is

class RunHeader(object):
    """
    Class representing the run header sub-block
    as specified in Corsika user manual, Table 7.

    """
    def __init__(self, subblock):
        self.id = subblock[0]
        self.run_number = subblock[1]
        self.date_start = subblock[2]
        self.version = subblock[3]

        self.observation_levels = subblock[4]
        self.observation_heights = numpy.array(subblock[5:15]) * units.cm

        self.spectral_slope = subblock[15]
        self.min_energy = subblock[16] * units.GeV
        self.max_energy = subblock[17] * units.GeV

        self.flag_EGS4 = subblock[18]
        self.flag_NKG = subblock[19]

        self.cutoff_hadrons = subblock[20] * units.GeV
        self.cutoff_muons = subblock[21] * units.GeV
        self.cutoff_electrons = subblock[22] * units.GeV
        self.cutoff_photons = subblock[23] * units.GeV

        self.C = numpy.array(subblock[24:74])

        self.x_inclined = subblock[74]
        self.y_inclined = subblock[75]
        self.z_inclined = subblock[76]
        self.theta_inclined = subblock[77]
        self.phi_inclined = subblock[78]

        self.n_showers = subblock[92]
        self.CKA = numpy.array(subblock[94:134])
        self.CETA = numpy.array(subblock[134:139])
        self.CSTRBA = numpy.array(subblock[139:150])

        self.x_scatter_Cherenkov = subblock[247]
        self.y_scatter_Cherenkov = subblock[248]
        self.atmospheric_layer_coundaries = numpy.array(subblock[249:254])
        self.a_atmospheric = numpy.array(subblock[254:259])
        self.b_atmospheric = numpy.array(subblock[259:264])
        self.c_atmospheric = numpy.array(subblock[264:269])
        self.NFLAIN = subblock[269]
        self.NFLDIF = subblock[270]
        self.NFLPIF = numpy.floor(subblock[271] / 100)
        self.NFLPI0 = subblock[271] % 100
        self.NFRAGM = numpy.floor(subblock[272] / 100)
        self.NFLCHE = subblock[272] % 100

    def thickness_to_height(self, thickness):
        """"Calculate height (cm) for given thickness (gramms/cm**2)"""

        if (thickness > 631.1):
            height = c_atmospheric[1] * numpy.log(b_atmospheric[1] /
                                                  (thickness - a_atmospheric[1]))
        elif (thickness > 271.7):
            height = c_atmospheric[2] * numpy.log(b_atmospheric[2] /
                                                  (thickness - a_atmospheric[2]))
        elif (thickness > 3.0395):
            height = c_atmospheric[3] * numpy.log(b_atmospheric[3] /
                                                  (thickness - a_atmospheric[3]))
        elif (thickness > 1.28292e-3):
            height = c_atmospheric[4] * numpy.log(b_atmospheric[4] /
                                                  (thickness - a_atmospheric[4]))
        else:
            height = (a_atmospheric[5] - thickness) / c_atmospheric[5]
        return height

    def height_to_thickness(self, height):
        """Thickness (gramms/cm**2) of atmosphere given a height (cm)"""

        if (height < 4.e5):
            thickness = (a_atmospheric[1] + b_atmospheric[1] *
                         numpy.exp(-height / c_atmospheric[1]))
        elif (height < 1.e6):
            thickness = (a_atmospheric[2] + b_atmospheric[2] *
                         numpy.exp(-height / c_atmospheric[2]))
        elif (height < 4.e6):
            thickness = (a_atmospheric[3] + b_atmospheric[3] *
                         numpy.exp(-height / c_atmospheric[3]))
        elif (height < 1.e7):
            thickness = (a_atmospheric[4] + b_atmospheric[4] *
                         numpy.exp(-height / c_atmospheric[4]))
        else:
            thickness = a_atmospheric[5] - height * c_atmospheric[5]

        return thickness

    def __str__(self):
        return textwrap.dedent("""\
            Run header:
                id: {run_n}
                date: {date}
                version: {version}
            """.format(run_n=self.run_number,
                       date=self.date_start,
                       version=self.version))


class EventHeader(object):
    """
    Class representing the event header sub-block
    as specified in Corsika user manual, Table 8.

    """
    def __init__(self, subblock):
        self.id = subblock[0]
        self.event_number = subblock[1]
        self.particle_id = subblock[2]
        self.particle = particles.id[subblock[2]]
        self.energy = subblock[3] * units.GeV
        self.starting_altitude = subblock[4] * units.g / units.cm2
        self.first_target = subblock[5]
        self.first_interaction_altitude = subblock[6] * units.cm2  # Bad units?
        self.p_x = subblock[7] * units.GeV
        self.p_y = subblock[8] * units.GeV
        self.p_z = - subblock[9] * units.GeV
        self.zenith = subblock[10] * units.rad
        self.azimuth = subblock[11] * units.rad

        self.n_seeds = subblock[12]
        self.seeds = numpy.array(zip(subblock[13:41:3],
                                     subblock[14:42:3],
                                     subblock[15:43:3]))

        self.run_number = subblock[43]
        self.date_start = subblock[44]
        self.version = subblock[45]

        self.n_observation_levels = subblock[46]
        self.observation_heights = numpy.array(subblock[47:57]) * units.cm

        self.spectral_slope = subblock[57]
        self.min_energy = subblock[58] * units.GeV
        self.max_energy = subblock[59] * units.GeV

        self.cutoff_hadrons = subblock[60] * units.GeV
        self.cutoff_muons = subblock[61] * units.GeV
        self.cutoff_electrons = subblock[62] * units.GeV
        self.cutoff_photons = subblock[63] * units.GeV

        self.NFLAIN = subblock[64]
        self.NFLDIF = subblock[65]
        self.NFLPI0 = subblock[66]
        self.NFLPIF = subblock[67]
        self.NFLCHE = subblock[68]
        self.NFRAGM = subblock[69]

        self.magnetic_field_x = subblock[70] * units.micro * units.tesla
        self.magnetic_field_z = subblock[71] * units.micro * units.tesla
        self.flag_EGS4 = subblock[72]
        self.flag_NKG = subblock[73]

        self.flag_hadron_model_low = subblock[74]
        self.flag_hadron_model_high = subblock[75]
        self.flag_Cherenkov = subblock[76]
        self.flag_neutrino = subblock[77]
        self.flag_curved = subblock[78]
        self.flag_computer = subblock[79]
        self.theta_min = subblock[80] * units.rad
        self.theta_max = subblock[81] * units.rad
        self.phi_min = subblock[82] * units.rad
        self.phi_max = subblock[83] * units.rad

        self.cherenkov_bunch = subblock[84]
        self.cherenkov_n_x = subblock[85]
        self.cherenkov_n_y = subblock[86]
        self.cherenkov_grid_x = subblock[87] * units.cm
        self.cherenkov_grid_y = subblock[88] * units.cm
        self.cherenkov_detector_x = subblock[89] * units.cm
        self.cherenkov_detector_y = subblock[90] * units.cm
        self.cherenkov_output_flag = subblock[91]

        self.array_rotation = subblock[92] * units.rad
        self.flag_extra_muon_information = subblock[93]

        self.multiple_scattering_step_length_factor = subblock[94]
        self.Cherenkov_wavelength_min = subblock[95] * units.nanometer
        self.Cherenkov_wavelength_max = subblock[96] * units.nanometer
        self.uses_of_Cherenkov_event = subblock[97]
        self.core_x = numpy.array(subblock[98:118]) * units.cm
        self.core_y = numpy.array(subblock[118:138]) * units.cm

        self.flag_SIBYLL = subblock[138]
        self.flag_SIBYLL_cross = subblock[139]
        self.flag_QGSJET = subblock[140]
        self.flag_QGSJET_cross = subblock[141]
        self.flag_DPMJET = subblock[142]
        self.flag_DPMJET_cross = subblock[143]
        self.flag_VENUS_cross = subblock[144]
        self.flag_muon_multiple = subblock[145]
        self.NKG_radial_range = subblock[146] * units.cm
        self.energy_fraction_thinning_hadronic = subblock[147]
        self.energy_fraction_thinning_em = subblock[148]
        self.weightlimit_thinning_hadronic = subblock[149]
        self.weightlimit_thinning_EM = subblock[150]
        self.radial_thinning_max_radius = subblock[151] * units.cm
        self.cone_inner_angle = subblock[152] * units.rad
        self.cone_outer_angle = subblock[153] * units.rad

        self.transistion_energy_high_low = subblock[154] * units.GeV
        self.flag_skimming_incidence = subblock[155]
        self.altitude_horizontal_shower_axis = subblock[156] * units.cm
        self.starting_height = subblock[157] * units.cm
        self.flag_charm = subblock[158]
        self.flag_hadron_origin = subblock[159]

        # CONEX
        self.min_vertical_depth_CONEX = subblock[160]
        self.high_threshold_hadron_CONEX = subblock[161]
        self.high_threshold_muon_CONEX = subblock[162]
        self.high_threshold_EM_CONEX = subblock[163]
        self.low_threshold_hadron_CONEX = subblock[164]
        self.low_threshold_muon_CONEX = subblock[165]
        self.low_threshold_EM_CONEX = subblock[166]
        self.flag_observation_level_curvature_CONEX = subblock[167]
        self.weightlimit_thinning_hadronic_CONEX = subblock[168]
        self.weightlimit_thinning_EM_CONEX = subblock[169]
        self.weightlimit_sampling_hadronic_CONEX = subblock[170]
        self.weightlimit_sampling_muons_CONEX = subblock[171]
        self.weightlimit_sampling_EM_CONEX = subblock[172]

    @property
    def hadron_model_low(self):
        hadron_models_low = {1: 'GHEISHA', 2: 'UrQMD', 3: 'FLUKA'}
        return hadron_models_low.get(self.flag_hadron_model_low, 'unknown')

    @property
    def hadron_model_high(self):
        hadron_models_high = {0: 'HDPM', 1: 'VENUS', 2: 'SIBYLL', 3: 'QGSJET',
                              4: 'DPMJET', 5: 'NEXUS', 6: 'EPOS'}
        return hadron_models_high.get(self.flag_hadron_model_high, 'unknown')

    @property
    def computer(self):
        computers = {3: 'UNIX', 4: 'Macintosh'}
        return computers.get(self.flag_computer, 'unknown')

    def __str__(self):
        return textwrap.dedent("""\
            Event header:
                id: {event_n}
                primary: {primary}
                energy: {energy} EeV
                direction: ({zenith}, {azimuth})
            """.format(event_n=self.event_number,
                       primary=self.particle_id,
                       energy=self.energy / units.EeV,
                       zenith=self.zenith / units.degree,
                       azimuth=self.azimuth / units.degree))


class RunEnd(object):
    """
    Class representing the run end sub-block
    as specified in Corsika user manual, Table 14.

    """
    def __init__(self, subblock):
        self.id = subblock[0]
        self.run_number = subblock[1]
        self.n_events_processed = subblock[2]

    def __str__(self):
        return textwrap.dedent("""\
            Run end:
                id: {run_number}
                events: {n_events}
            """.format(run_number=self.run_number,
                       n_events=self.n_events_processed))


class EventEnd(object):
    """
    Class representing the event end sub-block
    as specified in Corsika user manual, Table 13.

    """
    def __init__(self, subblock):
        self.id = subblock[0]
        self.event_number = subblock[1]

        self.n_photons_levels = subblock[2]
        self.n_electrons_levels = subblock[3]
        self.n_hadrons_levels = subblock[4]
        self.n_muons_levels = subblock[5]
        self.n_particles_levels = subblock[6]

        # NKG output
        self.NKG_lateral_1_x = numpy.array(subblock[7:28]) / units.cm2
        self.NKG_lateral_1_y = numpy.array(subblock[28:49]) / units.cm2
        self.NKG_lateral_1_xy = numpy.array(subblock[49:70]) / units.cm2
        self.NKG_lateral_1_yx = numpy.array(subblock[70:91]) / units.cm2

        self.NKG_lateral_2_x = numpy.array(subblock[91:112]) / units.cm2
        self.NKG_lateral_2_y = numpy.array(subblock[112:133]) / units.cm2
        self.NKG_lateral_2_xy = numpy.array(subblock[133:154]) / units.cm2
        self.NKG_lateral_2_yx = numpy.array(subblock[154:175]) / units.cm2

        self.NKG_electron_number = numpy.array(subblock[175:185])
        self.NKG_pseudo_age = numpy.array(subblock[185:195])
        self.NKG_electron_distances = numpy.array(subblock[195:205]) * units.cm
        self.NKG_local_pseudo_age_1 = numpy.array(subblock[205:215])

        self.NKG_level_height_mass = numpy.array(subblock[215:225])
        self.NKG_level_height_distance = numpy.array(subblock[225:235])
        self.NKG_distance_bins_local_pseudo_age = numpy.array(subblock[235:245]) * units.cm
        self.NKG_local_pseudo_age_2 = numpy.array(subblock[245:255])

        # Longitudinal distribution
        self.longitudinal_parameters = numpy.array(subblock[255:261])
        self.chi2_longitudinal_fit = subblock[261]

        self.n_photons_output = subblock[262]
        self.n_electrons_output = subblock[263]
        self.n_hadrons_output = subblock[264]
        self.n_muons_output = subblock[265]
        self.n_preshower_EM_particles = subblock[266]

    def __str__(self):
        return textwrap.dedent("""\
            Event end:
                id: {event_n}
                particles: {particles}
            """.format(event_n=self.event_number,
                       particles=self.n_particles_levels))


class ParticleData(object):
    """
    Class representing the particle data sub-block
    as specified in Corsika user manual, Table 10.

    """
    def __init__(self, subblock):
        self.description = subblock[0]
        self.p_x = subblock[1] * units.GeV
        self.p_y = subblock[2] * units.GeV
        self.p_z = - subblock[3] * units.GeV
        self.x = subblock[4] * units.cm
        self.y = subblock[5] * units.cm
        self.t = subblock[6] * units.ns  # or z for additional muon info

        self.id = numpy.floor(self.description / 1000)
        self.r = numpy.sqrt(self.x ** 2 + self.y ** 2)
        self.is_particle = 0 < self.id < 200
        self.particle = particles.id[self.id] if self.is_particle else None
        self.is_detectable = self.particle in ['positron', 'electron',
                                               'muon_p', 'muon_m']  # + 'gamma'

    @property
    def hadron_generation(self):
        return numpy.floor(self.description / 10) % 100

    @property
    def observation_level(self):
        return self.description % 10

    @property
    def phi(self):
        return numpy.arctan2(self.y, self.x)

    @property
    def is_nucleus(self):
        return 200 <= self.id < 9900 or self.id == 14

    @property
    def is_Cherenkov(self):
        return 9900 <= self.id

    @property
    def atomic_number(self):
        if self.is_nucleus:
            if self.id == 14:
                return 1
            else:
                return self.id % 100
        else:
            return -1

    @property
    def atom(self):
        if self.is_nucleus:
            return particles.atomic_number[self.atomic_number]
        else:
            return None

    def __str__(self):
        return textwrap.dedent("""\
            Particle:
                description: {description}
                momentum: {momentum} GeV
                position: {position} m
                time: {time} ns
            """.format(description=self.description,
                       momentum=(self.p_x / units.GeV,
                                 self.p_y / units.GeV,
                                 self.p_z / units.GeV),
                       position=(self.x / units.m, self.y / units.m),
                       time=self.t / units.ns))


class CherenkovData(object):
    """
    Class representing the cherenkov photon sub-block
    as specified in Corsika user manual, Table 11.

    The number of CherenkovData records in a sub-block depends on
    compilation options.

    """
    def __init__(self, subblock):
        self.photons_in_bunch = subblock[0]
        self.x = subblock[1] * units.cm
        self.y = subblock[2] * units.cm
        self.u = subblock[3]
        self.v = subblock[4]
        self.t = subblock[5] * units.ns
        self.production_height = subblock[6] * units.cm

    def __str__(self):
        return textwrap.dedent("""\
            Cherenkov:
                n: {n}
                position: {position} m
                (u,v): {direction}
                time: {time} ns
                prod. height: {height} m
            """.format(n=self.photons_in_bunch,
                       direction=(self.u, self.v),
                       position=(self.x / units.m, self.y / units.m),
                       time=self.t / units.ns,
                       height=self.production_height / units.m))


# THIN versions

class FormatThin(Format):
    """
    Class containing the format information of the thinned file
    as specified in Corsika user manual, Section 10.2.2.

    """
    def __init__(self):
        super(FormatThin, self).__init__()

        # one block contains 6552 fields plus one header and one end field
        self.block_format = '6554f'
        self.block_size = struct.calcsize(self.block_format)

        # Each block contains 21 sub-blocks
        # Each sub-block consists of 312 fields
        # the first of which _might_ be a string id.
        self.subblock_format = '4s311f'
        self.subblock_size = struct.calcsize(self.subblock_format)

        # Each particle record sub-block contains a fixed
        # number of particle records
        # With the thinned option, each of these is 8 fields long
        # for a total of 39 records per sub block
        self.particle_format = '8f'
        self.particle_size = struct.calcsize(self.particle_format)


class ParticleDataThin(ParticleData):
    """
    Class representing the thinned particle data sub-block
    as specified in Corsika user manual, Table 10.

    """
    def __init__(self, subblock):
        self.weight = subblock[7]
        super(ParticleDataThin, self).__init__(subblock)

    def __str__(self):
        return textwrap.dedent("""\
            Particle:
                id: {description}
                momentum: {momentum} GeV
                position: {position} m
                time: {time} ns
                weight: {weight}
            """.format(description=self.description,
                       momentum=(self.p_x / units.GeV,
                                 self.p_y / units.GeV,
                                 self.p_z / units.GeV),
                       position=(self.x / units.m, self.y / units.m),
                       time=self.t / units.ns,
                       weight=self.weight))


class CherenkovDataThin(CherenkovData):
    """
    Class representing the thinned cherenkov photon sub-block
    as specified in Corsika user manual, Table 11.

    The number of CherenkovData records in a sub-block depends on
    compilation options.

    """
    def __init__(self, subblock):
        self.weight = subblock[7]
        super(CherenkovDataThin, self).__init__(subblock)

    def __str__(self):
        return textwrap.dedent("""\
            Cherenkov:
                n: {n}
                position: {position} m
                (u,v): {direction}
                time: {time} ns
                prod. height: {height} m
                weight: {weight}
            """.format(n=self.photons_in_bunch,
                       direction=(self.u, self.v),
                       position=(self.x / units.m, self.y / units.m),
                       time=self.t / units.ns,
                       height=self.production_height / units.m,
                       weight=self.weight))
