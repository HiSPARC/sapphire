"""
Use this for particle identification for particles in CORSIKA

Particle codes as specified in CORSIKA user manual, Table 4.

This relates to the Description in Particles as::

    particle id x 1000 + hadron generation x 10 + no. of observation level

So to find e+/e-:

.. code-block:: python

    from sapphire.corsika import particles

    particle.id = 2
    if particles.name(particle.id) in ['positron', 'electron']:
        pass

    particle.id = 1206
    if particle.id > 200:
        print 'atom: %s' % particles.ATOMIC_NUMBER[particle.id % 100]


"""
import re


def name(particle_id):
    """Get the name for a CORSIKA particle code

    :param particle_id: code for the particle
    :return: name of the particle. In case of atoms the weight is added
             to the name.

    """
    try:
        return ID[particle_id]
    except KeyError:
        return ATOMIC_NUMBER[particle_id % 100] + str(int(particle_id / 100))


def particle_id(name):
    """Get the CORSIKA particle code for a partice name

    :param name: name of the particle/atom, for atoms the mass
                 (neutrons + protons) can be appended to the name.
    :return: CORSIKA code for the particle.
             For atoms the code is: A x 100 + Z

    """
    for pid, particle_name in ID.iteritems():
        if name == particle_name:
            return pid
    for z, atom_name in ATOMIC_NUMBER.iteritems():
        # Note; the mass number assumes no neutrons. To get a correct
        # weight append the weight to the name, e.g. helium4 or carbon14
        if name == atom_name:
            return z * 100 + z
    atom = re.match("^([a-z]+)(\d+)$", name)
    if atom is not None:
        return int(atom.group(2)) * 100 + (particle_id(atom.group(1)) % 100)


ID = {1: 'gamma',
      2: 'positron',
      3: 'electron',
      4: 'neutrino',  # No longer used?
      5: 'muon_p',
      6: 'muon_m',
      7: 'pion_0',
      8: 'pion_p',
      9: 'pion_m',
      10: 'Kaon_0_long',
      11: 'Kaon_p',
      12: 'Kaon_m',
      13: 'neutron',
      14: 'proton',
      15: 'anti_proton',
      16: 'Kaon_0_short',
      17: 'eta',
      18: 'Lambda',
      19: 'Sigma_p',
      20: 'Sigma_0',
      21: 'Sigma_m',
      22: 'Xi_0',
      23: 'Xi_m',
      24: 'Omega_m',
      25: 'anti_neutron',
      26: 'anti_Lambda',
      27: 'anti_Sigma_m',
      28: 'anti_Sigma_0',
      29: 'anti_Sigma_p',
      30: 'anti_Xi_0',
      31: 'anti_Xi_p',
      32: 'anti_Omega_p',
      50: 'omega',
      51: 'rho_0',
      52: 'rho_p',
      53: 'rho_m',
      54: 'Delta_pp',
      55: 'Delta_p',
      56: 'Delta_0',
      57: 'Delta_m',
      58: 'anti_Delta_mm',
      59: 'anti_Delta_m',
      60: 'anti_Delta_0',
      61: 'anti_Delta_p',
      62: 'Kaon_star_0',
      63: 'Kaon_star_p',
      64: 'Kaon_star_m',
      65: 'anti_Kaon_star_0',
      66: 'electron_neutrino',
      67: 'anti_electron_neutrino',
      68: 'muon_neutrino',
      69: 'anti_muon_neutrino',

      71: 'eta__2_gamma',
      72: 'eta__3_pion_0',
      73: 'eta__pion_p_pion_m_pion_0',
      74: 'eta__pion_p_pion_m_gamma',
      75: 'additional_muon_p',
      76: 'additional_muon_m',

      85: 'decay_start_muon_p',
      86: 'decay_start_muon_m',

      95: 'decay_end_muon_p',
      96: 'decay_end_muon_m',

      116: 'D_0',
      117: 'D_p',
      118: 'anti_D_m',
      119: 'anti_D_0',
      120: 'D_p_short',
      121: 'anti_D_m_short',
      122: 'eta_c',
      123: 'D_star_0',
      124: 'D_star_p',
      125: 'anti_D_star_m',
      126: 'anti_D_star_0',
      127: 'D_star_p_short',
      128: 'anti_D_star_m_short',

      130: 'j_psi',
      131: 'tau_p',
      132: 'tau_m',
      133: 'tau_neutrino',
      134: 'anti_tau_neutrino',

      137: 'Lambda_c_p',
      138: 'Xi_c_p',
      139: 'Xi_c_0',
      140: 'Sigma_c_pp',
      141: 'Sigma_c_',
      142: 'Sigma_c_0',
      143: 'Xi_c_prime_p',
      144: 'Xi_c_prime_0',
      145: 'Omega_c_0',

      149: 'anti_Lambda_c_m',
      150: 'anti_Xi_c_m',
      151: 'anti_Xi_c_0',
      152: 'anti_Sigma_c_mm',
      153: 'anti_Sigma_c_m',
      154: 'anti_Sigma_c_0',
      155: 'anti_Xi_c_prime_m',
      156: 'anti_Xi_c_prime_0',
      157: 'anti_Omega_c_0',

      161: 'Sigma_c_star_pp',
      162: 'Sigma_c_star_p',
      163: 'Sigma_c_star_0',

      171: 'anti_Sigma_c_star_mm',
      172: 'anti_Sigma_c_star_m',
      173: 'anti_Sigma_c_star_0',

      176: 'B_0',
      177: 'B_p',
      178: 'anti_B_m',
      179: 'anti_B_0',
      180: 'B_s_0',
      181: 'anti_B_s_0',
      182: 'B_c_p',
      183: 'anti_B_c_m',
      184: 'Lambda_b_0',
      185: 'Sigma_b_m',
      186: 'Sigma_b_p',
      187: 'Xi_b_0',
      188: 'Xi_b_m',
      189: 'Omega_b_m',
      190: 'anti_Lambda_b_0',
      191: 'anti_Sigma_b_p',
      192: 'anti_Sigma_b_m',
      193: 'anti_Xi_b_0',
      194: 'anti_Xi_b_p',
      195: 'anti_Omega_b_p',

      # A x 100 + Z
      101: 'hydrogen',
      201: 'deuteron',
      301: 'tritium',
      302: 'helium3',
      402: 'alpha',
      703: 'lithium',
      904: 'beryllium',
      1105: 'boron',
      1206: 'carbon',
      1407: 'nitrogen',
      1608: 'oxygen',
      2713: 'aluminium',
      2814: 'silicon',
      3216: 'sulfur',
      4020: 'calcium',
      5626: 'iron',
      5828: 'nickel',

      9900: 'cherenkov_photons'}


# Z numbers
ATOMIC_NUMBER = {1: 'hydrogen',
                 2: 'helium',
                 3: 'lithium',
                 4: 'beryllium',
                 5: 'boron',
                 6: 'carbon',
                 7: 'nitrogen',
                 8: 'oxygen',
                 9: 'fluorine',
                 10: 'neon',
                 11: 'sodium',
                 12: 'magnesium',
                 13: 'aluminium',
                 14: 'silicon',
                 15: 'phosphorus',
                 16: 'sulfur',
                 17: 'chlorine',
                 18: 'argon',
                 19: 'potassium',
                 20: 'calcium',
                 21: 'scandium',
                 22: 'titanium',
                 23: 'vanadium',
                 24: 'chromium',
                 25: 'manganese',
                 26: 'iron',
                 27: 'cobalt',
                 28: 'nickel',
                 29: 'copper',
                 30: 'zinc',
                 31: 'gallium',
                 32: 'germanium',
                 33: 'arsenic',
                 34: 'selenium',
                 35: 'bromine',
                 36: 'krypton',
                 37: 'rubidium',
                 38: 'strontium',
                 39: 'yttrium',
                 40: 'zirconium',
                 41: 'niobium',
                 42: 'molybdenum',
                 43: 'technetium',
                 44: 'ruthenium',
                 45: 'rhodium',
                 46: 'palladium',
                 47: 'silver',
                 48: 'cadmium',
                 49: 'indium',
                 50: 'tin',
                 51: 'antimony',
                 52: 'tellurium',
                 53: 'iodine',
                 54: 'xenon',
                 55: 'caesium',
                 56: 'barium',
                 57: 'lanthanum',
                 58: 'cerium',
                 59: 'praseodym.',
                 60: 'neodymium',
                 61: 'promethium',
                 62: 'samarium',
                 63: 'europium',
                 64: 'gadolinium',
                 65: 'terbium',
                 66: 'dysprosium',
                 67: 'holmium',
                 68: 'erbium',
                 69: 'thulium',
                 70: 'ytterbium',
                 71: 'lutetium',
                 72: 'hafnium',
                 73: 'tantalum',
                 74: 'tungsten',
                 75: 'rhenium',
                 76: 'osmium',
                 77: 'iridium',
                 78: 'platinum',
                 79: 'gold',
                 80: 'mercury',
                 81: 'thallium',
                 82: 'lead',
                 83: 'bismuth',
                 84: 'polonium',
                 85: 'astatine',
                 86: 'radon',
                 87: 'francium',
                 88: 'radium',
                 89: 'actinium',
                 90: 'thorium',
                 91: 'protactin.',
                 92: 'uranium',
                 93: 'neptunium',
                 94: 'plutonium',
                 95: 'americium',
                 96: 'curium',
                 97: 'berkelium',
                 98: 'californium',
                 99: 'einsteinium'}


# From the CORSIKA `corsikaread.cpp` program
#
#    naming conventions of corsika particles:
#     1   gamma           24   Omega           64   K* -
#     2   positron        25   anti neutron    65   anti K* 0
#     3   electron        26   anti Lambda     66   electron neutrino
#     4   neutrino        27   anti Sigma -    67   electron anti neutrino
#     5   muon +          28   anti Sigma 0    68   muon neutrino
#     6   muon -          29   anti Sigma +    69   muon anti neutrino
#     7   pion 0          30   anti Xi 0       71   eta-> 2*gam
#     8   pion +          31   anti Xi +       72   eta-> 3*pi0
#     9   pion -          32   anti Omega      73   eta-> pi+ + pi- + pi0
#    10   Kaon 0 long     50   omega           74   eta-> pi+ + pi- + gam
#    11   Kaon            51   rho 0           201   Deuteron
#    12   Kaon -          52   rho +           301   Tritium
#    13   neutron         53   rho -           402   alpha
#    14   proton          54   Delta ++       1206   Carbon
#    15   anti proton     55   Delta +        1407   Nitrogen
#    16   Kaon 0 short    56   Delta 0        1608   Oxygen
#    17   eta (71..74)    57   Delta -        2713   Aluminium
#    18   Lambda          58   anti Delta --  2814   Silicon
#    19   Sigma +         59   anti Delta -   3216   Sulfur
#    20   Sigma 0         60   anti Delta 0   5626   Iron
#    21   Sigma -         61   anti Delta +   9900   Cherenkov photons
#    22   Xi 0            62   K* 0
#    23   Xi -            63   K* +
#   116   D 0            131   tau +           150   anti Xi c -
#   117   D +            132   tau -           151   anti Xi c 0
#   118   anti D -       133   tau neutrino    152   anti Sigma c --
#   119   anti D 0       134   anti tau neutr  153   anti Sigma c -
#   120   D s +          137   Lambda c +      154   anti Sigma c 0
#   121   anti D s -     138   Xi c +          155   anti Xi c prime -
#   122   eta c          139   Xi c 0          156   anti Xi c prime 0
#   123   D*0            140   Sigma c ++      157   anti Omega c 0
#   124   D*+            141   Sigma c +       161   Sigma c * ++
#   125   anti D*-       142   Sigma c 0       162   Sigma c * +
#   126   anti D*0       143   Xi c prime +    163   Sigma c * 0
#   127   D* s +         144   Xi c prime 0    171   anti Sigma c * --
#   128   anti D* s -    145   Omega c 0       172   anti Sigma c * -
#   130   J/psi          149   anti Lambda c-  173   anti Sigma c * 0
