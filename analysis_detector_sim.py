import tables
from itertools import combinations
import re

from pylab import *
from scipy.optimize import curve_fit
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

def plot_random_timing_errors(data, tablename, D, dt1=None, dt3=None,
                              dt4=None, limit=None):
    NS, NF, NT = 0, 0, 0

    N = 10

    phis, thetas = [], []
    dphis, dthetas = [], []
    table = data.getNode('/analysis', tablename)
    for event in table[:limit]:
        if min(event['n1'], event['n3'], event['n4']) >= D:
            NT += 1
            theta, phi = reconstruct_angle(event)
            if not isnan(theta) and not isnan(phi):
                NS += 1
                phis.append(phi)
                thetas.append(theta)

                t1 = event['t1']
                t3 = event['t3']
                t4 = event['t4']

                for i in range(N):
                    if dt1:
                        event['t1'] = t1 + normal(scale=dt1)
                    if dt3:
                        event['t3'] = t3 + normal(scale=dt3)
                    if dt4:
                        event['t4'] = t4 + normal(scale=dt4)
                    theta2, phi2 = reconstruct_angle(event)
                    dphis.append((phi, phi - phi2))
                    dthetas.append((theta, theta - theta2))
            else:
                NF += 1

    phis = array(phis)
    thetas = array(thetas)
    dthetas = array(dthetas)
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

    figure()
    plot(rad2deg(dphis[:,0]), rad2deg(dphis[:,1]), '.', ms=1.)
    xlabel("Azimuthal angle (deg)")
    ylabel("Error in azimuthal angle (deg)")
    title("Timing errors introduce azimuthal angle errors (rand)")
    
    figure()
    plot(rad2deg(dthetas[:,0]), rad2deg(dthetas[:,1]), '.', ms=1.)
    xlabel("Zenith angle (deg)")
    ylabel("Error in zenith angle (deg)")
    title("Timing errors introduce zenith angle errors (rand)")

    return phis, dphis, thetas, dthetas

def plot_estimate_timing_errors(dtv, whichdt, N=10, thetacond=None):
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    phis = zeros(shape=(N, N, N))
    thetas = zeros(shape=(N, N, N))
    phis_err = zeros(shape=(N, N, N))
    thetas_err = zeros(shape=(N, N, N))
    dt = 300. / (N - 1)
    for i, t1 in enumerate(linspace(-100, 200, N)):
        #t1 += uniform(-.5 * dt, .5 * dt)
        for j, t3 in enumerate(linspace(-100, 200, N)):
            #t3 += uniform(-.5 * dt, .5 * dt)
            for k, t4 in enumerate(linspace(-100, 200, N)):
                #t4 += uniform(-.5 * dt, .5 * dt)

                theta, phi = reconstruct_angle_dt(t1 - t3, t1 - t4)
                if thetacond:
                    if not thetacond[0] < rad2deg(theta) < thetacond[1]:
                        phis[i,j,k] = nan
                        thetas[i,j,k] = nan
                        continue

                if not isnan(theta) and not isnan(phi):
                    phis[i,j,k] = phi
                    thetas[i,j,k] = theta
                    if whichdt == 'dt1':
                        phis_err[i,j,k] = dphi_dt1(phi1, phi2, phi, t1,
                                                   t3, t4)
                    elif whichdt == 'dt3':
                        phis_err[i,j,k] = dphi_dt3(phi1, phi2, phi, t1,
                                                   t3, t4)
                        thetas_err[i,j,k] = dtheta_dt3(theta, phi1, phi2,
                                                       phi, t1, t3, t4)
                    elif whichdt == 'dt4':
                        phis_err[i,j,k] = dphi_dt4(phi1, phi2, phi, t1,
                                                   t3, t4)
                    else:
                        raise Exception("Unsupport dt type: %s" % whichdt)
                else:
                    phis[i,j,k] = nan
                    thetas[i,j,k] = nan

    d = array([(u, v) for u, v in zip(phis.flatten(), phis_err.flatten())
               if not isnan(u)])

    # Fit a sin curve to data
    f = lambda x, a, b: a * sin(x + b)
    popt, pcov = curve_fit(f, d[:,0], dtv * d[:,1])
    print "Fitting phis_err from %s (%.1f): parameters: %s" % (whichdt,
                                                               dtv, popt)
    phifitp = popt

    figure()
    plot(rad2deg(d[:,0]), rad2deg(dtv * d[:,1]), '.', ms=1.)
    plot(rad2deg(d[:,0]), rad2deg(-dtv * d[:,1]), '.', ms=1.)
    x = linspace(-pi, pi, 50)
    plot(rad2deg(x), rad2deg(f(x, popt[0], popt[1])))
    plot(rad2deg(x), -rad2deg(f(x, popt[0], popt[1])))

    xlabel("Azimuthal angle (deg)")
    ylabel("Error in azimuthal angle (deg)")
    title("Azimuthal angle error due to timing errors\n"
          "(%s, %.1f) (calc, %.1f) (theta, %s)" % (whichdt, dtv,
                                                   rad2deg(popt[0]),
                                                   thetacond))

    figure()
    hist(rad2deg(d[:,0]), bins=linspace(-180, 180, 90), histtype='step')

    if whichdt == 'dt3':
        d = array([(u, v, w) for u, v, w in zip(thetas.flatten(),
                                                thetas_err.flatten(),
                                                phis.flatten())
                   if not isnan(u)])
        figure()
        plot(rad2deg(d[:,0]), rad2deg(dtv * d[:,1]), '.', ms=1.)
        plot(rad2deg(d[:,0]), rad2deg(-dtv * d[:,1]), '.', ms=1.)
        xlabel("Zenith angle (deg)")
        ylabel("Error in zenith angle (deg)")
        title("Zenith angle error due to timing errors\n"
              "(%s, %.1f) (calc) (theta, %s)" % (whichdt, dtv, thetacond))

        figure()
        plot(rad2deg(d[:,2]), rad2deg(dtv * d[:,1]), '.', ms=1.)
        plot(rad2deg(d[:,2]), rad2deg(-dtv * d[:,1]), '.', ms=1.)
        xlabel("Azimuthal angle (deg)")
        ylabel("Error in zenith angle (deg)")
        title("Zenith angle error due to timing errors\n"
              "(%s, %.1f) (calc) (theta, %s)" % (whichdt, dtv, thetacond))

    return phis, phis_err, thetas, thetas_err, phifitp

def dphi_dt1(phi1, phi2, phi, t1, t3, t4):
    r1 = r2 = 10.

    return 1 / (1 + tan(phi) ** 2) * \
           (r2 * cos(phi2) - r1 * cos(phi1) + \
            (r2 * sin(phi2) - r1 * sin(phi1)) * tan(phi)) / \
           (r2 * (t3 - t1) * sin(phi2) - r1 * (t4 - t1) * sin(phi1))

def dphi_dt3(phi1, phi2, phi, t1, t3, t4):
    r1 = r2 = 10.

    return 1 / (1 + tan(phi) ** 2) * \
           (-r2 * (sin(phi2) * tan(phi) + cos(phi2))) / \
           (r2 * (t3 - t1) * sin(phi2) - r1 * (t4 - t1) * sin(phi1))

def dphi_dt4(phi1, phi2, phi, t1, t3, t4):
    r1 = r2 = 10.

    return 1 / (1 + tan(phi) ** 2) * \
           (r1 * (cos(phi1) + sin(phi1) * tan(phi))) / \
           (r2 * (t3 - t1) * sin(phi2) - r1 * (t4 - t1) * sin(phi1))

def dtheta_dt3(theta, phi1, phi2, phi, t1, t3, t4):
    r1 = 10.
    # Since times are in ns, c must be in m / ns
    c = 3.00e8 / 1e9

    return 1 / sqrt(1 - sin(theta) ** 2) * \
           (c + r1 * sin(theta) * sin(phi - phi1) * \
            dphi_dt3(phi1, phi2, phi, t1, t3, t4)) / \
           (r1 * cos(phi - phi1))

def plot_total_phi_error(poptdt1, poptdt3, poptdt4):
    x = linspace(-pi, pi, 100)
    y1 = abs(poptdt1[0] * sin(x + poptdt1[1]))
    y3 = abs(poptdt3[0] * sin(x + poptdt3[1]))
    y4 = abs(poptdt4[0] * sin(x + poptdt4[1]))
    Y = sqrt(y1 ** 2 + y3 ** 2 + y4 ** 2)
    figure()
    plot(rad2deg(x), rad2deg(y1), label="Error from dt1")
    plot(rad2deg(x), rad2deg(y3), label="Error from dt3")
    plot(rad2deg(x), rad2deg(y4), label="Error from dt4")
    plot(rad2deg(x), rad2deg(Y), label="Total error")
    xlabel("Azimuthal angle (deg)")
    ylabel("Error in azimuthal angle (deg)")
    title("Total error in azimuthal angle")
    legend(loc="lower right")

def mytest(N=10, phicond=None):
    phi1 = calc_phi(1, 3)
    phi2 = calc_phi(1, 4)

    phis = zeros(shape=(N, N, N))
    thetas = zeros(shape=(N, N, N))
    phis_err = zeros(shape=(N, N, N))
    dt = 300. / (N - 1)
    for i, t1 in enumerate(linspace(-100, 200, N)):
        t1 += uniform(-.5 * dt, .5 * dt)
        for j, t3 in enumerate(linspace(-100, 200, N)):
            t3 += uniform(-.5 * dt, .5 * dt)
            for k, t4 in enumerate(linspace(-100, 200, N)):
                t4 += uniform(-.5 * dt, .5 * dt)

                theta, phi = reconstruct_angle_dt(t1 - t3, t1 - t4)

                if not isnan(theta) and not isnan(phi):
                    phis[i,j,k] = phi
                    thetas[i,j,k] = theta
                    phis_err[i,j,k] = dphi_dt4(phi1, phi2, phi, t1, t3,
                                               t4)
                    if phicond:
                        if phicond[0] < rad2deg(phi) < phicond[1]:
                            print phi, phis_err[i,j,k], (t3 - t1) / (t4 - t1)
                else:
                    phis[i,j,k] = nan
                    thetas[i,j,k] = nan
    return phis, phis_err



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

    N = 100
    #phis, dphis, thetas, dthetas = \
    #    plot_random_timing_errors(data, 'angle_23', 2, dt3=1.3,
    #                              limit=10000)
    #phis, phis_err, thetas, thetas_err, poptdt1 = \
    #    plot_estimate_timing_errors(1.3, 'dt1', N=N)#, thetacond=(18, 28))
    #phis, phis_err, thetas, thetas_err, poptdt3 = \
    #    plot_estimate_timing_errors(1.3, 'dt3', N=N)#, thetacond=(18, 28))
    #phis, phis_err, thetas, thetas_err, poptdt4 = \
    #    plot_estimate_timing_errors(1.3, 'dt4', N=N)#, thetacond=(18, 28))
    #plot_total_phi_error(poptdt1, poptdt3, poptdt4)
