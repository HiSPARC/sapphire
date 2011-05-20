import tables

try:
    data
except NameError:
    data = tables.openFile('data-e15.h5', 'r')

try:
    R, T
except NameError:
    events = data.root.simulations.angle_0

    r_list, t_list = [], []
    for event in events:
        if not event['id'] % 10:
            r0 = event['r']
            fi0 = event['phi']
        elif event['pid'] in [-2, 2, -3, 3]:
            r1 = event['r']
            fi1 = event['phi']

            # uit mail jos
            #x=r0*cos(fi0)+r1*cos(fi1)
            #y=r0*sin(fi0)+r1*sin(fi1)
            #r=sqrt(x*x+y*y)
            ###
            # correcte versie (zie pamflet)
            r = r1

            r_list.append(r)
            t_list.append(event['time'])
    R = array(r_list)
    T = array(t_list)

figure()

r = linspace(0, 150, 50)
mr, tm, dtmin, dtmax = [], [], [], []
for r0, r1 in zip(r[:-1], r[1:]):
    t = compress((r0 <= R) & (R < r1), T)
    tm.append(mean(t))
    dtmin.append(sqrt(mean(compress(t - t.mean() < 0, t - t.mean()) ** 2)))
    dtmax.append(sqrt(mean(compress(t - t.mean() > 0, t - t.mean()) ** 2)))
    mr.append((r0 + r1) / 2.)
#errorbar(mr, tm, yerr=(dtmin, dtmax), drawstyle='steps-mid', capsize=0,
#         label="simulation (csv)")
errorbar(mr, tm, drawstyle='steps-mid', capsize=0,
         label="simulation (csv)")
print "CSV max:", max(r_list)

l = data.root.showers.zenith0.leptons
mr, tm, dtmin, dtmax = [], [], [], []
for r0, r1 in zip(r[:-1], r[1:]):
    t = l.readWhere('(r0 <= core_distance) & (core_distance < r1)')['arrival_time']
    tm.append(mean(t))
    dtmin.append(sqrt(mean(compress(t - t.mean() < 0, t - t.mean()) ** 2)))
    dtmax.append(sqrt(mean(compress(t - t.mean() > 0, t - t.mean()) ** 2)))
    mr.append((r0 + r1) / 2.)
#errorbar(mr, tm, yerr=(dtmin, dtmax), drawstyle='steps-mid', capsize=0,
#         label="shower")
errorbar(mr, tm, drawstyle='steps-mid', capsize=0,
         label="shower")

l = data.root.analysis.angle_0
mr, tm, dtmin, dtmax = [], [], [], []
for r0, r1 in zip(r[:-1], r[1:]):
    t = l.readWhere('(r0 <= r) & (r < r1)')['t1']
    t = compress(-isnan(t), t)
    tm.append(mean(t))
    dtmin.append(sqrt(mean(compress(t - t.mean() < 0, t - t.mean()) ** 2)))
    dtmax.append(sqrt(mean(compress(t - t.mean() > 0, t - t.mean()) ** 2)))
    mr.append((r0 + r1) / 2.)
#errorbar(mr, tm, yerr=(dtmin, dtmax), drawstyle='steps-mid', capsize=0,
#         label="simulation (scint1 FP))")
errorbar(mr, tm, drawstyle='steps-mid', capsize=0,
         label="simulation (scint1 FP))")

xlabel("Core distance (m)")
ylabel("Mean arrival time (ns)")
title(r"$\theta = 0^\circ, E = 1\,\mathrm{PeV}$")
legend()
