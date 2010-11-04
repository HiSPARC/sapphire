import tables
from itertools import combinations
import re

from pylab import *
from tikz_plot import tikz_2dhist


DETECTORS = [(0., 5.77, 'UD'), (0., 0., 'UD'), (-5., -2.89, 'LR'),
             (5., -2.89, 'LR')]


def plot_all_ring_timings(data):
    """Plot timing histograms for various core distances"""

    plot_ring_timings(data, [(0, 4), (4, 20), (20, 40), (40, 80),
                             (80, 120)], normed=False, binstep=.1)
    plot_ring_timings(data, [(40, 50), (50, 60), (60, 70), (70, 80)],
                      normed=True, binstep=.5)

def plot_ring_timings(data, rings, normed, binstep):
    """Plot timing histograms for various core distances"""

    figure()
    for r0, r1 in rings:
        t = []
        events = data.root.analysis.readWhere('(r0 <= r) & (r < r1)')
        times = events['times'].T
        for s1, s2 in combinations(range(4), 2):
            t.extend(times[s1] - times[s2])
        t = [x for x in t if not isnan(x)]
        hist(t, bins=arange(-20, 20, binstep), histtype='step',
             normed=normed, label="%.1f < r < %.1f" % (r0, r1))
    legend()
    title("Time differences between scintillator events")
    xlabel("time (ns)")
    ylabel("count")

def plot_reconstructed_angles(data, tablename, THETA, D, binning=False,
                              randomize_binning=False, N=None):
    """Reconstruct angles from simulation for minimum particle density"""

    match = re.search('_size([0-9]+)', tablename)
    if match:
        R = int(match.group(1))
    else:
        R = 10

    NS = 0
    NF = 0
    NT = 0
    alphas = []
    thetas = []
    phis = []
    opening_angles = []
    table = data.getNode('/analysis', tablename)
    for event in table[:N]:
        if min(event['n1'], event['n3'], event['n4']) >= D:
            NT += 1
            # Do we need to bin timing data?
            if binning is not False:
                event['t1'] = floor(event['t1'] / binning) * binning 
                event['t2'] = floor(event['t2'] / binning) * binning 
                event['t3'] = floor(event['t3'] / binning) * binning 
                event['t4'] = floor(event['t4'] / binning) * binning 
                # Do we need to randomize inside a bin?
                if randomize_binning is True:
                    event['t1'] += uniform(0, binning)
                    event['t2'] += uniform(0, binning)
                    event['t3'] += uniform(0, binning)
                    event['t4'] += uniform(0, binning)

            theta, phi = reconstruct_angle(event, R)
            if not isnan(theta) and not isnan(phi):
                NS += 1
                alpha = event['alpha']
                alphas.append(alpha)
                thetas.append(theta)
                phis.append(phi)
                opening_angle = arccos(sin(theta) * sin(THETA) *
                                       cos(phi - -alpha) + cos(theta) *
                                       cos(THETA))
                opening_angles.append(opening_angle)
            else:
                NF += 1

    opening_angles = array(opening_angles)
    alphas = array(alphas)
    thetas = array(thetas)
    phis = array(phis)

    subtitle = "theta: %.1f degrees, size: %d meters, binsize: %r, randomized: %r" % (rad2deg(THETA), R, binning, randomize_binning)

    figure()
    plot(-rad2deg(alphas), rad2deg(phis), '.', ms=1.)
    axis('tight')
    title("Reconstructed azimuthal angle from simulation (D >= %d)\n%s" % (D, subtitle))
    xlabel("Simulated angle (degrees)")
    ylabel("Reconstructed angle (degrees)")
    savefig('plots/auto-azimuth-TH%d-D%d-R%d-b%d-rb%s.pdf' % (rad2deg(THETA), D, R, binning, randomize_binning))

    fn = 'plots/auto-azimuth-TH%d-D%d-R%d-b%d-rb%s.tikz' % (rad2deg(THETA), D, R, binning, randomize_binning)
    #tikz_2dhist(fn, -rad2deg(alphas), rad2deg(phis), bins=(72,72),
    tikz_2dhist(fn, -rad2deg(alphas), rad2deg(phis), bins=(9,9),
                use_log=False)

    figure()
    hist(rad2deg(thetas), bins=linspace(0, 90, 200), histtype='step')
    title("Reconstructed zenith angle from simulation (D >= %d)\n%s" % (D, subtitle))
    xlabel("Reconstructed angle (degrees)")
    ylabel("Count")
    savefig('plots/auto-zenith-TH%d-D%d-R%d-b%d-rb%s.pdf' % (rad2deg(THETA), D, R, binning, randomize_binning))

    figure()
    n, bins, patches = hist(rad2deg(opening_angles),
                            bins=linspace(0, 120, 200), histtype='step')
    res = get_resolution(n, bins)
    axvspan(xmin=0, xmax=res, color='blue', alpha=.2)
    title("Opening angle of reconstructed and simulated angles (D >= %d)\n%s" % (D, subtitle))
    xlabel("Opening angle (degrees)")
    ylabel("Count")
    figtext(.65, .8, "resolution: %.2f deg\nfailed: %5.1f %%" %
            (res, 100. * NF / NT))
    savefig('plots/auto-opening-TH%d-D%d-R%d-b%d-rb%s.pdf' % (rad2deg(THETA), D, R, binning, randomize_binning))

    print
    print
    print "Angle reconstruction (D >= %d)" % D
    print "Size of station: %d meters" % R
    print "Simulated zenith angle: %.1f degrees" % rad2deg(THETA)
    if binning is not False:
        print "Binning of timings was used with binsize: %f ns" % binning
        if randomize_binning is True:
            print "Timing values were randomized inside a bin."
    else:
        print "Unbinned timings used."
    print "Total of %d (%d) events" % (len(thetas),
                                       len(table))
    print "Total number reconstructions:        %6d" % NT
    print "Number of succesful reconstructions: %6d (%5.1f %%)" % \
        (NS, 100. * NS / NT)
    print "Number of failed reconstructions:    %6d (%5.1f %%)" % \
        (NF, 100. * NF / NT)
    print "Angle resolution (68%% integrated): %.2f degrees" % res
    print

def reconstruct_angle(event, R=10):
    """Reconstruct angles from a single event"""

    dt1 = event['t1'] - event['t3']
    dt2 = event['t1'] - event['t4']

    return reconstruct_angle_dt(dt1, dt2, R)

def reconstruct_angle_dt(dt1, dt2, R=10):
    """Reconstruct angle given time differences"""

    c = 3.00e+8

    r1 = R
    r2 = R 

    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                  (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
    theta = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
    theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))

    return theta, phi

def calc_phi(s1, s2):
    x1, y1 = DETECTORS[s1 - 1][:2]
    x2, y2 = DETECTORS[s2 - 1][:2]

    return arctan2((y2 - y1), (x2 - x1))

def get_resolution(n, bins):
    """Get angle resolution from histogram values

    Resolution is defined as the opening angle which contains 68 % of all
    events.

    """
    total = sum(n)
    N = 0 
    for c, res in zip(n, bins[1:]):
        N += c
        if 1. * N / total >= .68:
            break
    return res

def plot_random_angles(N):
    DT = 10 / 3e8 * 1e9
    ts1 = uniform(-DT / 2, DT / 2, N)
    ts3 = uniform(-DT / 2, DT / 2, N)
    ts4 = uniform(-DT / 2, DT / 2, N)

    thetas, phis = [], []
    for t1, t3, t4 in zip(ts1, ts3, ts4):
        t1 = floor(t1 / 2.5) * 2.5
        t3 = floor(t3 / 2.5) * 2.5
        t4 = floor(t4 / 2.5) * 2.5
        event = dict(t1=t1, t3=t3, t4=t4)
        theta, phi = reconstruct_angle(event)
        thetas.append(theta)
        phis.append(phi)

    print "NaNs in phi: %d (of %d)" % (len([x for x in phis if isnan(x)]),
                                       len(phis))
    nnan = len([x for x in thetas if isnan(x)])
    print "NaNs in theta: %d (of %d) (%3.2f)" % (nnan, len(thetas),
                                                 1. * nnan / len(thetas) * 100)

    figure()
    hist([x for x in phis if not isnan(x)], bins=200, histtype='step')
    figure()
    hist([x for x in thetas if not isnan(x)], bins=200, histtype='step')

    return thetas, phis

def plot_all_reconstructed_angles(data):
    """Generate plots used in pamflet"""

    # zenith 0 degrees
    plot_reconstructed_angles(data, 'angle_0', 0, D=2)

    # zenith 5 degrees, D=1,2,4
    kwargs = dict(data=data, tablename='angle_5', THETA=deg2rad(5))
    plot_reconstructed_angles(D=1, **kwargs)
    plot_reconstructed_angles(D=2, **kwargs)
    plot_reconstructed_angles(D=4, **kwargs)

    # zenith 22.5 degrees, D=1,2,4
    kwargs = dict(data=data, tablename='angle_23', THETA=pi / 8)
    plot_reconstructed_angles(D=1, **kwargs)
    plot_reconstructed_angles(D=2, **kwargs)
    plot_reconstructed_angles(D=4, **kwargs)

    # zenith 35 degrees, D=1,2,4
    kwargs = dict(data=data, tablename='angle_35', THETA=deg2rad(35))
    plot_reconstructed_angles(D=1, **kwargs)
    plot_reconstructed_angles(D=2, **kwargs)
    plot_reconstructed_angles(D=4, **kwargs)

    # SPECIALS
    # zenith 22.5, D=2, sizes=5,20
    kwargs = dict(data=data, THETA=pi / 8, D=2)
    plot_reconstructed_angles(tablename='angle_23_size5', **kwargs)
    plot_reconstructed_angles(tablename='angle_23_size20', **kwargs)

    # zenith 22.5, D=2, binnings
    kwargs = dict(data=data, tablename='angle_23', THETA=pi / 8, D=2)
    plot_reconstructed_angles(binning=1, randomize_binning=True, **kwargs)
    plot_reconstructed_angles(binning=2.5, randomize_binning=False, **kwargs)
    plot_reconstructed_angles(binning=2.5, randomize_binning=True, **kwargs)
    plot_reconstructed_angles(binning=5, randomize_binning=True, **kwargs)

def plot_estimate_timing_errors(data, tablename, D, N, limit=None):
    NS, NF, NT = 0, 0, 0

    phis = []
    dphis = []
    table = data.getNode('/analysis', tablename)
    for event in table[:limit]:
        if min(event['n1'], event['n3'], event['n4']) >= D:
            NT += 1
            theta, phi = reconstruct_angle(event)
            if not isnan(theta) and not isnan(phi):
                NS += 1
                phis.append(phi)

                t1 = event['t1']
                t3 = event['t3']
                t4 = event['t4']

                for i in range(N):
                    dt1, dt3, dt4 = normal(scale=1.3, size=3)
                    event['t1'] = t1 + dt1
                    event['t3'] = t3 + dt3
                    event['t4'] = t4 + dt4
                    theta2, phi2 = reconstruct_angle(event)
                    dphis.append((phi, phi - phi2))
            else:
                NF += 1

    phis = array(phis)
    # Make sure all dphis are within [-pi, pi)
    dphis = (array(dphis) + pi) % (2 * pi) - pi

    print "Total of %d (%d) events" % (len(phis),
                                       len(table))
    print "Total number reconstructions:        %6d" % NT
    print "Number of succesful reconstructions: %6d (%5.1f %%)" % \
        (NS, 100. * NS / NT)
    print "Number of failed reconstructions:    %6d (%5.1f %%)" % \
        (NF, 100. * NF / NT)
    print

    plot(dphis[:,0], dphis[:,1], '.', ms=1.)
    xlabel("Azimuthal angle (rad)")
    ylabel("Error in azimuthal angle (rad)")
    title("Timing errors introduce azimuthal angle errors")

    return phis, dphis


if __name__ == '__main__':
    # invalid values in arcsin will be ignored (nan handles the situation
    # quite well)
    np.seterr(invalid='ignore')

    try:
        data
    except NameError:
        data = tables.openFile('data-e15.h5', 'r')

    # For pamflet:
    #plot_all_reconstructed_angles(data)

    phis, dphis = plot_estimate_timing_errors(data, 'angle_23', 2, 10,
                                              limit=10000)
