from sapphire.analysis.direction_reconstruction import DirectionReconstruction
from sapphire.clusters import ScienceParkCluster

TIMING = 5.5

rec = DirectionReconstruction
cluster = ScienceParkCluster([501, 503, 506])
r1, phi1 = cluster.calc_r_and_phi_for_stations(0, 1)
r2, phi2 = cluster.calc_r_and_phi_for_stations(0, 2)
phis = linspace(-pi, pi, 50)

dt = sqrt(mean(rec.rel_phi_errorsq(pi / 8, phis, phi1, phi2, r1, r2)))
dphi = dt * TIMING
print rad2deg(dphi)

dt = sqrt(mean(rec.rel_theta1_errorsq(pi / 8, phis, phi1, phi2, r1, r2)))
dtheta = dt * TIMING
print rad2deg(dtheta)

print rad2deg(sqrt((dphi * sin(pi / 8)) ** 2 + dtheta ** 2))
