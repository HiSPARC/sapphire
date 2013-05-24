"""
Use this for particle identification for particles in corsika

Particle codes as specified in Corsika user manual, Table 4.

This relates to the Description in Particles as:
particle id x 1000 + hadron generation x 10 + no. of observation level

So to find e+/e-:

    if floor(particle.fDescription / 1000) in [positron, electron]:
        pass

"""

gamma = 1
positron = 2
electron = 3
neutrino = 4  # No longer used?
muon_p = 5
muon_m = 6
pion_0 = 7
pion_p = 8
pion_m = 9
kaon_0_long = 10
kaon_p = 11
kaon = kaon_p
kaon_m = 12
neutron = 13
proton = 14
anti_proton = 15
kaon_0_short = 16
eta = 17
lambda_ = 18
sigma_p = 19
sigma_0 = 20
sigma_m = 21
xi_0 = 22
xi_m = 23
omega = 24
anti_neutron = 25
anti_lambda = 26
anti_sigma_m = 27
anti_sigma_0 = 28
anti_sigma_p = 29
anti_xi_0 = 30
anti_xi_p = 31
anti_omega = 32
w = 50
rho_0 = 51
rho_p = 52
rho_m = 53
delta_pp = 54
delta_p = 55
delta_0 = 56
delta_m = 57
anti_delta_mm = 58
anti_delta_m = 59
anti_delta_0 = 60
anti_delta_p = 61
k_star_0 = 62
k_star_p = 63
k_star_m = 64
k_star_0_anti = 65
electron_neutrino =  66
electron_anti_neutrino =  67
muon_neutrino =  68
muon_anti_neutrino =  69
eta__2_gamma =  71
eta__3_pion_0 =  72
eta__pion_p_pion_m_pion_0 =  73
eta__pion_p_pion_m_gamma =  74

D_0 = 116
D_p = 117
anti_D_m = 118
anti_D_0 = 119
D_p_short = 120
anti_D_m_short = 121
eta_c = 122
D_star_0 = 123
D_star_p = 124
anti_D_star_m = 125
anti_D_star_0 = 126
D_star_p_short = 127
anti_D_star_m_short= 128
j_psi = 130
tau_p = 131
tau_m = 132
tau_neutrino = 133
anti_tau_neutrino = 134
Lambda_c_p = 137
Xi_c_p = 138
Xi_c_0 = 139
Sigma_c_pp = 140
Sigma_c_ = 141
Sigma_c_0 = 142
Xi_c_prime_p = 143
Xi_c_prime_0 = 144
Omega_c_0 = 145
anti_Lambda_c_m = 149
anti_Xi_c_m = 150
anti_Xi_c_0 = 151
anti_Sigma_c_mm = 152
anti_Sigma_c_m = 153
anti_Sigma_c_0 = 154
anti_Xi_c_prime_m = 155
anti_Xi_c_prime_0 = 156
anti_Omega_c_0 = 157
Sigma_c_star_pp = 161
Sigma_c_star_p = 162
Sigma_c_star_0 = 163
anti_Sigma_c_star_mm = 171
anti_Sigma_c_star_m = 172
anti_Sigma_c_star_0 = 173

# A x 100 + Z
deuteron = 201
tritium = 301
alpha = 402
carbon = 1206
nitrogen = 1407
oxygen = 1608
aluminium = 2713
silicon = 2814
sulfur = 3216
iron = 5626

cherenkov_photons = 9900


# From the CORSIKA `corsikaread.cpp` program
#
#    naming conventions of corsika particles:
#     1   gamma           24   omega           64   k* -
#     2   positron        25   anti neutron    65   k* 0 anti
#     3   electron        26   anti lambda     66   electron neutrino
#     4   neutrino        27   anti sigma -    67   electron anti neutrino
#     5   muon +          28   anti sigma 0    68   muon neutrino
#     6   muon -          29   anti sigma +    69   muon anti neutrino
#     7   pion 0          30   anti xi 0       71   eta-> 2*gam
#     8   pion +          31   anti xi +       72   eta-> 3*pi0
#     9   pion -          32   anti omega      73   eta-> pi+ + pi- + pi0
#    10   kaon 0 long     50   w               74   eta-> pi+ + pi- + gam
#    11   kaon            51   rho 0           201   deuteron
#    12   kaon -          52   rho +           301   tritium
#    13   neutron         53   rho -           402   alpha
#    14   proton          54   delta ++       1206   carbon
#    15   anti proton     55   delta +        1407   nitrogen
#    16   kaon 0 short    56   delta 0        1608   oxygen
#    17   eta (71..74)    57   delta -        2713   aluminium
#    18   lambda          58   anti delta --  2814   silicon
#    19   sigma +         59   anti delta -   3216   sulfur
#    20   sigma 0         60   anti delta 0   5626   iron
#    21   sigma -         61   anti delta +   9900   cherenkov photons
#    22   xi 0            62   k* 0
#    23   xi -            63   k* +
#   116   D 0            131   tau +           150   anti Xi c -
#   117   D +            132   tau -           151   anti Xi c 0
#   118   anti D -       133   tau neutrino    152   anti Sigma c --
#   119   anti D 0       134   anti tau neutr  153   anti Sigma c -
#   120   D+ short       137   Lambda c+       154   anti Sigma c 0
#   121   anti D- short  138   Xi c +          155   anti Xi c prime -
#   122   eta c          139   Xi c 0          156   anti Xi c prime 0
#   123   D*0            140   Sigma c ++      157   anti Omega c 0
#   124   D*+            141   Sigma c +       161   Sigma c * ++
#   125   anti D*-       142   Sigma c 0       162   Sigma c * +
#   126   anti D*0       143   Xi c prime +    163   Sigma c * 0
#   127   D*+ short      144   Xi c prime 0    171   anti Sigma c * --
#   128   anti D*- short 145   Omega c 0       172   anti Sigma c * -
#   130   j/psi          149   anti Lambda c-  173   anti Sigma c * 0
