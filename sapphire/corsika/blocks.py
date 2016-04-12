"""
Classes corresponding to CORSIKA blocks and sub-blocks

The classes in this module correspond one-to-one with the sub-blocks
(and blocks) as specified in the CORSIKA users manual.

Author: Javier Gonzalez <jgonzalez@ik.fzk.de>
"""

import struct
import math

import numpy

import units
import particles
try:
    from numba import jit
except ImportError:
    def jit(func):
        return func


# All sizes are in bytes

class Format(object):

    """The binary format information of the file.

    As specified in the CORSIKA user manual, Section 10.2.1.

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
        self.fields_per_particle = 7
        self.particle_format = '%df' % self.fields_per_particle
        self.particle_size = struct.calcsize(self.particle_format)
        self.particles_per_subblock = 39

        # Full particle sub block
        self.particles_format = (self.particle_format *
                                 self.particles_per_subblock)
        self.particles_size = self.particle_size * self.particles_per_subblock


# From here on, things should not depend on the field size as everything is

class RunHeader(object):

    """The run header sub-block

    As specified in the CORSIKA user manual, Table 7.

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
        self.atmospheric_layer_boundaries = numpy.array(subblock[249:254])
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
        """"Calculate height (m) for given thickness (gramms/cm**2)

        As specified in the CORSIKA user manual, Appendix D.

        """
        a, b, c = self.a_atmospheric, self.b_atmospheric, self.c_atmospheric

        layers = [l * units.cm for l in self.atmospheric_layer_boundaries]

        if thickness > self.height_to_thickness(layers[1]):
            height = c[0] * math.log(b[0] / (thickness - a[0]))
        elif thickness > self.height_to_thickness(layers[2]):
            height = c[1] * math.log(b[1] / (thickness - a[1]))
        elif thickness > self.height_to_thickness(layers[3]):
            height = c[2] * math.log(b[2] / (thickness - a[2]))
        elif thickness > self.height_to_thickness(layers[4]):
            height = c[3] * math.log(b[3] / (thickness - a[3]))
        else:
            height = (a[4] - thickness) * c[4] / b[4]

        return height * units.cm

    def height_to_thickness(self, height):
        """Thickness (gramms/cm**2) of atmosphere given a height (m)

        As specified in the CORSIKA user manual, Appendix D.

        """
        height = height * units.m / units.cm
        a, b, c = self.a_atmospheric, self.b_atmospheric, self.c_atmospheric

        if height < self.atmospheric_layer_boundaries[1]:
            # 0-4 km
            thickness = a[0] + b[0] * math.exp(-height / c[0])
        elif height < self.atmospheric_layer_boundaries[2]:
            # 4-10 km
            thickness = a[1] + b[1] * math.exp(-height / c[1])
        elif height < self.atmospheric_layer_boundaries[3]:
            # 10-40 km
            thickness = a[2] + b[2] * math.exp(-height / c[2])
        elif height < self.atmospheric_layer_boundaries[4]:
            # 40-100 km
            thickness = a[3] + b[3] * math.exp(-height / c[3])
        else:
            # >100 km
            thickness = a[4] - b[4] * height / c[4]

        return thickness


class EventHeader(object):

    """The event header sub-block

    As specified in the CORSIKA user manual, Table 8.

    """

    def __init__(self, subblock):
        self.id = subblock[0]
        self.event_number = subblock[1]
        self.particle_id = subblock[2]
        self.particle = particles.name(subblock[2])
        self.energy = subblock[3] * units.GeV
        self.starting_altitude = subblock[4] * units.g / units.cm2
        self.first_target = subblock[5]
        self.first_interaction_altitude = subblock[6] * units.cm
        self.p_x = subblock[7] * units.GeV
        self.p_y = subblock[8] * units.GeV
        self.p_z = - subblock[9] * units.GeV  # Same direction as axis
        self.zenith = subblock[10] * units.rad

        # CORSIKA coordinate conventions are shown in Figure 1 of the manual.
        # CORSIKA defines azimuth as the direction the shower points to,
        # HiSPARC defines azimuth as the direction the shower comes from.
        # CORSIKA allows azimuths in [-2pi, 2pi], HiSPARC uses [-pi, pi).
        # CORSIKA defines North as 0 rad, HiSPARC defines East as 0 rad.
        # So finally we need to subtract pi/2 rad from the azimuth and
        # normalize its range.
        azimuth_corsika = subblock[11] * units.rad
        azimuth = azimuth_corsika - (math.pi / 2.)
        if azimuth >= math.pi:
            self.azimuth = azimuth - (2 * math.pi)
        elif azimuth < -math.pi:
            self.azimuth = azimuth + (2 * math.pi)
        else:
            self.azimuth = azimuth

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


class RunEnd(object):

    """The run end sub-block

    As specified in the CORSIKA user manual, Table 14.

    """

    def __init__(self, subblock):
        self.id = subblock[0]
        self.run_number = subblock[1]
        self.n_events_processed = subblock[2]


class EventEnd(object):

    """The event end sub-block

    As specified in the CORSIKA user manual, Table 13.

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
        self.NKG_distance_bins_local_pseudo_age = \
            numpy.array(subblock[235:245]) * units.cm
        self.NKG_local_pseudo_age_2 = numpy.array(subblock[245:255])

        # Longitudinal distribution
        self.longitudinal_parameters = numpy.array(subblock[255:261])
        self.chi2_longitudinal_fit = subblock[261]

        self.n_photons_output = subblock[262]
        self.n_electrons_output = subblock[263]
        self.n_hadrons_output = subblock[264]
        self.n_muons_output = subblock[265]
        self.n_preshower_EM_particles = subblock[266]


@jit
def particle_data(subblock):
    """Get particle data.

    As specified in the CORSIKA user manual, Table 10. High-performing
    version of the ParticleData class, but without all the easy-to-use
    attribute access.

    Transformations are needed for the x, y and p_z values from CORSIKA.
    CORSIKA coordinate conventions are mentioned in Figure 1 and Table 10.

    :return: tuple with p_x, p_y, p_z, x, y, t, id, r, hadron_generation,
             observation_level, phi data.

    """
    # These three are subject to coordinate transformations
    x_corsika = subblock[4] * units.cm
    y_corsika = subblock[5] * units.cm
    p_z_corsika = subblock[3] * units.GeV

    description = int(subblock[0])
    p_x = subblock[1] * units.GeV
    p_y = subblock[2] * units.GeV
    p_z = -p_z_corsika
    x = -y_corsika
    y = x_corsika
    t = subblock[6] * units.ns  # or z for additional muon info

    id = description / 1000
    hadron_generation = description / 10 % 100
    observation_level = description % 10

    r = math.sqrt(x ** 2 + y ** 2)
    phi = math.atan2(y, x)

    return (p_x, p_y, p_z, x, y, t, id, r, hadron_generation,
            observation_level, phi)


class ParticleData(object):

    """The particle data sub-block

    As specified in the CORSIKA user manual, Table 10.

    """

    def __init__(self, subblock):
        self.p_x, self.p_y, self.p_z, self.x, self.y, self.t, self.id, \
            self.r, self.hadron_generation, self.observation_level, \
            self.phi = particle_data(subblock)

    @property
    def is_detectable(self):
        """Get True or False if the particle is detectable

        Note: gamma particles are currently not included.

        """
        return self.particle in ['positron', 'electron', 'muon_p', 'muon_m']

    @property
    def particle(self):
        return particles.name(self.id) if self.is_particle else None

    @property
    def is_particle(self):
        return 0 < self.id < 200

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
            return particles.name(self.atomic_number)
        else:
            return None


class CherenkovData(object):

    """The cherenkov photon sub-block

    As specified in CORSIKA user manual, Table 11.

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


# THIN versions

class FormatThin(Format):

    """The format information of the thinned file

    As specified in CORSIKA user manual, Section 10.2.2.

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

    """The thinned particle data sub-block

    As specified in the CORSIKA user manual, Table 10.

    """

    def __init__(self, subblock):
        self.weight = subblock[7]
        super(ParticleDataThin, self).__init__(subblock)


class CherenkovDataThin(CherenkovData):

    """The thinned cherenkov photon sub-block

    As specified in CORSIKA user manual, Table 11.

    The number of CherenkovData records in a sub-block depends on
    compilation options.

    """

    def __init__(self, subblock):
        self.weight = subblock[7]
        super(CherenkovDataThin, self).__init__(subblock)
