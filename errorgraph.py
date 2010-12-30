import tables
data = tables.openFile('data-e15.h5', 'r')

figure()
l = data.root.showers.zenith0.leptons
r = linspace(0, 150, 50)
mr, tm, dtmin, dtmax = [], [], [], []
for r0, r1 in zip(r[:-1], r[1:]):
    t = l.readWhere('(r0 <= core_distance) & (core_distance < r1)')['arrival_time']
    tm.append(mean(t))
    dtmin.append(sqrt(mean(compress(t - t.mean() < 0, t - t.mean()) ** 2)))
    dtmax.append(sqrt(mean(compress(t - t.mean() > 0, t - t.mean()) ** 2)))
    mr.append((r0 + r1) / 2.)
    
errorbar(mr, tm, yerr=(dtmin, dtmax), drawstyle='steps-mid', capsize=0)
xlabel("Core distance (m)")
ylabel("Mean arrival time (ns)")
title(r"$\theta = 0^\circ, E = 1\,\mathrm{PeV}$")

#figure()
l = data.root.analysis.angle_0
r = linspace(0, 150, 50)
mr, tm, dtmin, dtmax = [], [], [], []
for r0, r1 in zip(r[:-1], r[1:]):
    t = l.readWhere('(r0 <= r) & (r < r1)')['t1']
    t = compress(-isnan(t), t)
    tm.append(mean(t))
    dtmin.append(sqrt(mean(compress(t - t.mean() < 0, t - t.mean()) ** 2)))
    dtmax.append(sqrt(mean(compress(t - t.mean() > 0, t - t.mean()) ** 2)))
    mr.append((r0 + r1) / 2.)
    
errorbar(mr, tm, yerr=(dtmin, dtmax), drawstyle='steps-mid', capsize=0)
xlabel("Core distance (m)")
ylabel("Mean arrival time (ns)")
title(r"$\theta = 0^\circ, E = 1\,\mathrm{PeV}$")
