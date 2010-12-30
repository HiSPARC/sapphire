import tables

try:
    data
except NameError:
    data = tables.openFile('data-e15.h5', 'r')

figure()
r = linspace(0, 100, 50)
events = data.root.analysis.angle_0
mr, tm, dt = [], [], []
for r0, r1 in zip(r[:-1], r[1:]):
    sel = events.readWhere('(r0 <= r) & (r < r1)')
    t1 = sel['t1'] - sel['t2']
    t3 = sel['t3'] - sel['t2']
    t4 = sel['t4'] - sel['t2']
    t = array([t1.flatten(), t3.flatten(), t4.flatten()]).flatten()
    t = compress(-isnan(t), t)
    tm.append(mean(t))
    dt.append(std(t))
    mr.append((r0 + r1) / 2.)
errorbar(mr, tm, yerr=(dtmin, dtmax), drawstyle='steps-mid', capsize=0,
         label="simulation")

xlabel("Core distance (m)")
ylabel("Mean arrival time (ns)")
title(r"$\theta = 0^\circ, E = 1\,\mathrm{PeV}$")
legend()

figure()
plot(mr, dt, 'o')
xlabel("Core distance (m)")
ylabel("Std dev arrival times (ns)")
title(r"$\theta = 0^\circ, E = 1\,\mathrm{PeV}$")
