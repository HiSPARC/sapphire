from __future__ import division

import tables
from scipy.special import gamma

data = tables.openFile('kascade.h5', 'r')
events = data.root.efficiency.events

def hist_shower_size():
    figure()
    hist(log10(events[:]['k_Num_e']), bins=100, histtype='step')
    xlabel("$\mathrm{lg}N_e$")
    ylabel("Count")
    title("Shower sizes in KASCADE data")
    axvline(4.6)
    annotate('Claimed KASCADE threshold', (4.6, 25000), (6, 25000),
             arrowprops={'arrowstyle': '->'})

def plot_ldf_from_data(size=(4.0, 4.2), limit=None):
    q = '(%f <= log10(k_Num_e)) & (log10(k_Num_e) < %f)' % size
#    figure()
#    sel = events.readWhere(q)
#    print len(sel)
#    plot(sel[:limit]['core_dist'],
#         sel['k_dens_e'][:limit,1] * cos(sel[:limit]['k_theta']), ',', label="All")
#    sel = sel.compress(sel[:]['k_theta'] < deg2rad(5))
#    plot(sel[:limit]['core_dist'],
#         sel['k_dens_e'][:limit,1] * cos(sel[:limit]['k_theta']), ',',
#         label=r"$\theta \leq 5^\circ$")
#    x = logspace(-1, 3)
#    plot(x, nkg(x, 10 ** size[0], .9))
#    plot(x, nkg(x, 10 ** size[1], .9))
#    xscale('log')
#    yscale('log')
#    xlabel("Core distance [m]")
#    ylabel("Electron density [$m^{-2}$]")
#    title(r"$%.1f \leq \mathrm{lg}(N_e) < %.1f$" % size)
#    legend(loc='best')

    figure()
    sel = events.readWhere(q)
    plot(sel[:]['core_dist'] * cos(sel[:]['k_theta']),
         sel['k_dens_e'][:,1] * cos(sel[:]['k_theta']), ',', label="All")
    sel = sel.compress(sel[:]['k_theta'] < deg2rad(5))
    plot(sel[:]['core_dist'] * cos(sel[:]['k_theta']),
         sel['k_dens_e'][:,1] * cos(sel[:]['k_theta']), ',',
         label=r"$\theta \leq 5^\circ$")
    x = logspace(-1, 3)
    plot(x, nkg(x, 10 ** size[0], .9))
    plot(x, nkg(x, 10 ** size[1], .9))
    xscale('log')
    yscale('log')
    xlabel(r"Core distance ($\cos\theta$-corrected) [m]")
    ylabel("Electron density [$m^{-2}$]")
    title(r"$%.1f \leq \mathrm{lg}(N_e) < %.1f$" % size)
    legend(loc='best')

def nkg(r, Ne, s):
    a = 1.5
    b = 3.6
    r0 = 40
    return Ne * c(s, a, b, r0) * (r / r0) ** (s - a) * \
           (1 + r / r0) ** (s - b)

def c(s, a, b, r0):
    return gamma(b - s) / (2 * pi * r0 ** 2 * gamma(s - a + 2) * \
                           gamma(a + b - 2 * s - 2))


if __name__ == '__main__':
    hist_shower_size()
    #plot_ldf_from_data()
    #plot_ldf_from_data((4.9, 5.1))
