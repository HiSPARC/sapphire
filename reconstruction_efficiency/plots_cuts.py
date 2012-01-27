from __future__ import division

import tables
from scipy.optimize import curve_fit, fmin_powell, fmin

from nkg import nkg
from build_dataset import DETECTORS
from plots import plot_hist_ratio


def main():
    #plot_ldfs()
    #fit_ldfs_core_pos()
    #hist_densities()
    #trigger_levels()
    #cuts_orientation_reconstruction()
    poisson_Ne()
    #poisson_s()
    pass

def trig_poisson(x):
    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    return 1 - (pd0 + 4 * pd1)          # Trigger probability

def plot_ldfs():
    events = data.root.efficiency.events

    q = "(k_Num_e > 10 ** 4.6) & (k_Num_e < 10 ** 7) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)

    q += "& (k_Num_e >= 10 ** 4.99) & (k_Num_e <= 10 ** 5.01)"

    sel = events.readWhere(q)

    clf()
    loglog(sel['core_dist'] / cos(sel['k_zenith']), sel['k_dens_e'][:,1],
           ',')

    x = sel['core_dist'] / cos(sel['k_zenith'])
    y = sel['k_dens_e'][:,1]
    Ne = 10 ** 5
    f = lambda x, a: nkg(x, Ne, a)
    popt, pcov = curve_fit(f, x, y, p0=(1.,))
    s, = popt
    r = logspace(-1, 3)
    plot(r, nkg(r, Ne, s), label="lg Ne = %f, s = %f" % (log10(Ne), s))

    xlabel("Core distance [m]")
    ylabel("Electron density [m$^{-2}$]")
    title("Lateral distribution")
    legend(loc='best')
    savefig('plots/ldfs_cuts.pdf')

def fit_ldfs_core_pos():
    events = data.root.efficiency.events

    q = "(k_Num_e > 10 ** 4.6) & (k_Num_e < 10 ** 7) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)

    #q += "& (k_Num_e >= 10 ** 4.99) & (k_Num_e <= 10 ** 5.01)"

    sel = events.readWhere(q)


    x = sel['k_core_pos'][:,0]
    y = sel['k_core_pos'][:,1]
    dens = sel['k_dens_e'][:,1]
    Ne = 10 ** 5
    s = 1.04
    f = lambda (a, b): sum((dens - nkg(sqrt((x - a) ** 2 + (y - b) ** 2),
                                       Ne, s))
                           ** 2)
    print f((65, 20.82))
    def g(x):
        print x, f(x)
    #sx, sy = fmin(f, (50, 50), xtol=1e-7, ftol=1e-7, callback=g)
    #print sx, sy

    figure()
    px = linspace(-100, 100, 100)
    py = linspace(-100, 100, 100)
    pz = array([[f((u, v)) for u in px] for v in py])
    contourf(px, py, pz, 30)
    dx = [u[0] for u in DETECTORS]
    dy = [u[1] for u in DETECTORS]
    scatter(dx, dy, c='r')
    savefig('plots/fit_core_pos.pdf')

    figure()
    px = linspace(59, 71, 100)
    py = linspace(14, 25, 100)
    pz = array([[f((u, v)) for u in px] for v in py])
    contourf(px, py, pz, 30)
    dx = [u[0] for u in DETECTORS]
    dy = [u[1] for u in DETECTORS]
    scatter(dx, dy, c='r')
    savefig('plots/fit_core_pos-zoom.pdf')

def hist_densities():
    events = data.root.efficiency.events

    q = "(k_Num_e > 10 ** 4.6) & (k_Num_e < 10 ** 7) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)

    #q += "& (k_Num_e >= 10 ** 4.99) & (k_Num_e <= 10 ** 5.01)"

    sel = events.readWhere(q)
    tsel = events.readWhere(q + '& (self_triggered == True)')

    clf()
    ax1 = subplot(111)
    ax1.hist(sel['k_cosdens_charged'][:,1], histtype='step',
             bins=linspace(0, 10, 100), normed=True)

    x = linspace(0, 10)
    # Poisson probabilities
    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit
    ptrig = 1 - (pd0 + 4 * pd1)         # Trigger probability

    ax1.plot(x, 1 - pd0, label="1 or more hit")
    ax1.plot(x, 1 - (pd0 + 4 * pd1), label="2 or more hit")
    ax1.plot(x, 1 - (pd0 + 4 * pd1 + 6 * pd2), label="3 or more hit")
    ax1.plot(x, 1 - (pd0 + 4 * pd1 + 6 * pd2 + 4 * pd3), label="4 hit")

    bins = linspace(0, 10, 51)
    plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
                        bins=bins, label="Data")
    legend()

    ax2 = ax1.twinx()
    x = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
    y = []
    for d0, d1 in zip(bins[:-1], bins[1:]):
        d = sel['k_cosdens_charged'][:,1]
        dsel = sel.compress((d0 <= d) & (d < d1))
        y.append(median(mean(dsel['pulseheights'], 1)))
    ax2.plot(x, y, label="Median pulseheight")

    savefig("plots/tries-poissons.pdf")

def trigger_levels():
    events = data.root.efficiency.events

    q = "(k_Num_e > 10 ** 4.6) & (k_Num_e < 10 ** 7) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)

    sel = events.readWhere(q)

    clf()
    bins = linspace(0, 1, 51)
    x = bins
    # Poisson probabilities
    p0 = exp(-.5 * x)                   # zero particles
    pp = 1 - p0                         # NOT zero particles
    pd0 = p0 ** 4 * pp ** 0             # 0 detectors hit
    pd1 = p0 ** 3 * pp ** 1             # 1 detector hit
    pd2 = p0 ** 2 * pp ** 2             # 2 detectors hit 
    pd3 = p0 ** 1 * pp ** 3             # 3 detectors hit
    pd4 = p0 ** 0 * pp ** 4             # 4 detectors hit
    ptrig = 1 - (pd0 + 4 * pd1)         # Trigger probability

    plot(x, ptrig, label="Poisson trigger")

    for level in arange(0, 300, 20):
        level = level / .57
        tsel = sel.compress(sum(sel['pulseheights'] >= level, 1) >= 2)
        plot_hist_ratio(tsel['k_cosdens_charged'],
                        sel['k_cosdens_charged'], bins=bins,
                        label='%s' % level)
    legend()
    savefig("plots/trigger-levels.pdf")

def cuts_orientation_reconstruction():
    global sel
    events = data.root.efficiency.events

    q = "(k_Num_e > 10 ** 4.6) & (k_Num_e < 10 ** 7) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)

    sel = events.readWhere(q)
    print len(sel)
    sel = sel.compress((sel['k_dens_e'][:,1] <= 0.2))
    print len(sel)
    #sel = sel.compress((-isnan(sel['h_phi'])))
    #print len(sel)

    clf()
    plot(sel['k_phi'], sel['h_phi'], ',')
    xlabel("$\phi_K$ [$^\degree$]")
    ylabel("$\phi_H$ [$^\degree$]")

    #plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
    #                    bins=bins, label="Data")

    return

    H, xedges, yedges = histogram2d(rad2deg(sel['k_phi']),
                                    rad2deg(sel['h_phi']), bins=200)
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2
    gca().set_axis_bgcolor('b')
    contourf(x, y, H.T, 20)
    colorbar()
    xlabel(r"$\phi_K\,(^\circ)$")
    ylabel(r"$\phi_H\,(^\circ)$")

def poisson_Ne():
    events = data.root.efficiency.events

    bins = linspace(0, 10)

    clf()

    q = "(k_Num_e >= 10 ** 4.6) & (k_Num_e < 10 ** 7) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)
    sel = events.readWhere(q)
    tsel = sel.compress(sel['self_triggered'] == True)
    plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
                        bins=bins, label="4.6 <= lg Ne < 7")

    q = "(k_Num_e < 10 ** 4.6) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)
    sel = events.readWhere(q)
    tsel = sel.compress(sel['self_triggered'] == True)
    plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
                        bins=bins, label="lg Ne < 4.6")

    x = bins
    plot(x, trig_poisson(x), label="Poisson")
    xlabel("Charged particle density [^m{-2}]")
    ylabel("Trigger probability")
    legend()
    savefig('plots/poisson_Ne.eps')

    axvline(0.2)

def poisson_s():
    events = data.root.efficiency.events

    bins = linspace(0, 10)

    clf()

    q = "(k_Num_e >= 10 ** 4.6) & (k_Num_e < 10 ** 7) & " \
        "(core_dist_center <= 90) & (k_theta <= %f)" % deg2rad(30)

    sel = events.readWhere(q)
    sel = sel.compress(sel['s'] <= .9)
    tsel = sel.compress(sel['self_triggered'] == True)
    plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
                        bins=bins, label="s <= 0.9")

    sel = events.readWhere(q)
    sel = sel.compress(sel['s'] <= 1)
    tsel = sel.compress(sel['self_triggered'] == True)
    plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
                        bins=bins, label="s <= 1")

    sel = events.readWhere(q)
    sel = sel.compress((1 <= sel['s']) & (sel['s'] < 1.5))
    tsel = sel.compress(sel['self_triggered'] == True)
    plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
                        bins=bins, label="1 <= s < 1.5")

    sel = events.readWhere(q)
    sel = sel.compress((sel['s'] > 1.5) & (sel['s'] < 2.))
    tsel = sel.compress(sel['self_triggered'] == True)
    plot_hist_ratio(tsel['k_cosdens_charged'], sel['k_cosdens_charged'],
                        bins=bins, label="s > 1.5")
    print len(sel)

    x = bins
    plot(x, trig_poisson(x), label="Poisson")
    xlabel("Charged particle density [^m{-2}]")
    ylabel("Trigger probability")
    legend()
    savefig('plots/poisson_s.eps')


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'r')

    main()
