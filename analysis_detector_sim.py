import tables
from itertools import combinations

from pylab import *


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

def plot_reconstructed_angles(data, D, N=None):
    """Reconstruct angles from simulation for minimum particle density"""

    THETA = pi / 8

    alphas = []
    thetas = []
    phis = []
    opening_angles = []
    for event in data.root.analysis[:N]:
        if min(event['n1'], event['n3'], event['n4']) >= D:
            #event['t1'] = floor(event['t1'] / 2.5) * 2.5 + uniform(0, 2.5)
            #event['t2'] = floor(event['t2'] / 2.5) * 2.5 + uniform(0, 2.5)
            #event['t3'] = floor(event['t3'] / 2.5) * 2.5 + uniform(0, 2.5)
            #event['t4'] = floor(event['t4'] / 2.5) * 2.5 + uniform(0, 2.5)
            theta, phi = reconstruct_angle(event)
            if not isnan(theta) and not isnan(phi):
                alpha = event['alpha']
                alphas.append(alpha)
                thetas.append(theta)
                phis.append(phi)
                opening_angle = arccos(sin(theta) * sin(THETA) *
                                       cos(phi - -alpha) + cos(theta) *
                                       cos(THETA))
                opening_angles.append(opening_angle)

    figure()
    plot([-x for x in alphas], phis, '.', ms=1.)
    axis('tight')
    title("Reconstructed azimuthal angle from simulation (D >= %d)" % D)
    xlabel("Simulated angle (radians)")
    ylabel("Reconstructed angle (radians)")
    figure()
    hist(thetas, bins=linspace(0, pi/2, 200), histtype='step')
    title("Reconstructed zenith angle from simulation (D >= %d)" % D)
    xlabel("Reconstructed angle (radians)")
    ylabel("Count")
    figure()
    hist(opening_angles, bins=200, histtype='step')
    title("Opening angle of reconstructed and simulated shower angle")
    xlabel("Opening angle (radians)")
    ylabel("Count")

    print "Angle reconstruction (D >= %d)" % D
    print "Total of %d (%d) events" % (len(thetas),
                                       len(data.root.analysis))
    opening_angles = [rad2deg(x) for x in opening_angles]
    print "Opening angle mean:    %.2f degrees" % mean(opening_angles)
    print "Opening angle std dev: %.2f degrees" % std(opening_angles)

def reconstruct_angle(event):
    """Reconstruct angles from a single event"""

    c = 3.00e+8

    r1 = 10
    r2 = 10

    dt1 = event['t1'] - event['t3']
    dt2 = event['t1'] - event['t4']

    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                  (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
    theta = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))
    theta2 = arcsin(c * dt2 * 1e-9 / (r2 * cos(phi - phi2)))
    #print 80 * '-'
    #if abs(theta - theta2) > 1e-5:
    #    print 80 * '!'
    #    print theta, theta2
    #if r1 < 0:
    #    print "R1", r1
    #if cos(phi - phi1) < 0:
    #    print "cos", cos(phi - phi1)
    #if dt1 < 0:
    #    print "dt1", dt1
    #if theta < 0:
    #    print "theta:", theta

    return theta, phi

def calc_phi(s1, s2):
    x1, y1 = DETECTORS[s1 - 1][:2]
    x2, y2 = DETECTORS[s2 - 1][:2]

    return arctan2((y2 - y1), (x2 - x1))

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


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('simulation-e15-angled.h5', 'r')

    #plot_all_ring_timings(data)
    _ip.magic("time plot_reconstructed_angles(data, D=1)")
    #_ip.magic("time plot_reconstructed_angles(data, D=4)")
    #thetas, phis = plot_random_angles(100000)
