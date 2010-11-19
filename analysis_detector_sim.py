import tables
from itertools import combinations
import re

from pylab import *
from scipy.optimize import curve_fit
from tikz_plot import tikz_2dhist


DETECTORS = [(0., 5.77, 'UD'), (0., 0., 'UD'), (-5., -2.89, 'LR'),
             (5., -2.89, 'LR')]


class ReconstructedEvent(tables.IsDescription):
    r = tables.Float32Col()
    phi = tables.Float32Col()
    alpha = tables.Float32Col()
    t1 = tables.Float32Col()
    t2 = tables.Float32Col()
    t3 = tables.Float32Col()
    t4 = tables.Float32Col()
    n1 = tables.UInt16Col()
    n2 = tables.UInt16Col()
    n3 = tables.UInt16Col()
    n4 = tables.UInt16Col()
    sim_theta = tables.Float32Col()
    sim_phi = tables.Float32Col()
    r_theta = tables.Float32Col()
    r_phi = tables.Float32Col()
    D = tables.UInt16Col()
    size = tables.UInt8Col()
    bin = tables.Float32Col()
    bin_r = tables.BoolCol()


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
    """Calculate angle between detectors (phi1, phi2)"""

    x1, y1 = DETECTORS[s1 - 1][:2]
    x2, y2 = DETECTORS[s2 - 1][:2]

    return arctan2((y2 - y1), (x2 - x1))

def dphi_dt1(phi, theta, phi1, phi2, r1=10, r2=10):
    r1 = r2 = 10.

    return 1 / (1 + tan(phi) ** 2) * \
           (r2 * cos(phi2) - r1 * cos(phi1) + \
            (r2 * sin(phi2) - r1 * sin(phi1)) * tan(phi)) / \
           (r2 * (t3 - t1) * sin(phi2) - r1 * (t4 - t1) * sin(phi1))

def do_full_reconstruction(data, tablename):
    """Do a reconstruction of all simulated data and store results"""

    try:
        data.createGroup('/', 'reconstructions', "Angle reconstructions")
    except tables.NodeError:
        pass

    try:
        data.removeNode('/reconstructions', tablename)
    except tables.NoSuchNodeError:
        pass

    table = data.createTable('/reconstructions', tablename,
                             ReconstructedEvent, "Reconstruction data")

    # zenith 0 degrees
    kwargs = dict(data=data, tablename='angle_0', THETA=0, dest=table)
    reconstruct_angles(D=1, **kwargs)
    reconstruct_angles(D=2, **kwargs)
    reconstruct_angles(D=3, **kwargs)
    reconstruct_angles(D=4, **kwargs)
    reconstruct_angles(D=5, **kwargs)

    # zenith 5 degrees, D=1,2,3,4,5
    kwargs = dict(data=data, tablename='angle_5', THETA=deg2rad(5),
                  dest=table)
    reconstruct_angles(D=1, **kwargs)
    reconstruct_angles(D=2, **kwargs)
    reconstruct_angles(D=3, **kwargs)
    reconstruct_angles(D=4, **kwargs)
    reconstruct_angles(D=5, **kwargs)

    # zenith 22.5 degrees, D=1,2,3,4,5
    kwargs = dict(data=data, tablename='angle_23', THETA=pi / 8,
                  dest=table)
    reconstruct_angles(D=1, **kwargs)
    reconstruct_angles(D=2, **kwargs)
    reconstruct_angles(D=3, **kwargs)
    reconstruct_angles(D=4, **kwargs)
    reconstruct_angles(D=5, **kwargs)

    # zenith 35 degrees, D=1,2,3,4,5
    kwargs = dict(data=data, tablename='angle_35', THETA=deg2rad(35),
                  dest=table)
    reconstruct_angles(D=1, **kwargs)
    reconstruct_angles(D=2, **kwargs)
    reconstruct_angles(D=3, **kwargs)
    reconstruct_angles(D=4, **kwargs)
    reconstruct_angles(D=5, **kwargs)

    # SPECIALS
    # zenith 22.5, sizes=5,20, D=1,2,3,4,5
    kwargs = dict(data=data, THETA=pi / 8, dest=table)
    reconstruct_angles(tablename='angle_23_size5', D=1, **kwargs)
    reconstruct_angles(tablename='angle_23_size5', D=2, **kwargs)
    reconstruct_angles(tablename='angle_23_size5', D=3, **kwargs)
    reconstruct_angles(tablename='angle_23_size5', D=4, **kwargs)
    reconstruct_angles(tablename='angle_23_size5', D=5, **kwargs)
    reconstruct_angles(tablename='angle_23_size20', D=1, **kwargs)
    reconstruct_angles(tablename='angle_23_size20', D=2, **kwargs)
    reconstruct_angles(tablename='angle_23_size20', D=3, **kwargs)
    reconstruct_angles(tablename='angle_23_size20', D=4, **kwargs)
    reconstruct_angles(tablename='angle_23_size20', D=5, **kwargs)

    # zenith 22.5, binnings, D=1,2,3,4,5
    kwargs = dict(data=data, tablename='angle_23', THETA=pi / 8,
                  dest=table)
    kwargs['randomize_binning'] = False
    reconstruct_angles(binning=1, D=1, **kwargs)
    reconstruct_angles(binning=1, D=2, **kwargs)
    reconstruct_angles(binning=1, D=3, **kwargs)
    reconstruct_angles(binning=1, D=4, **kwargs)
    reconstruct_angles(binning=1, D=5, **kwargs)
    reconstruct_angles(binning=2.5, D=1, **kwargs)
    reconstruct_angles(binning=2.5, D=2, **kwargs)
    reconstruct_angles(binning=2.5, D=3, **kwargs)
    reconstruct_angles(binning=2.5, D=4, **kwargs)
    reconstruct_angles(binning=2.5, D=5, **kwargs)
    reconstruct_angles(binning=5, D=1, **kwargs)
    reconstruct_angles(binning=5, D=2, **kwargs)
    reconstruct_angles(binning=5, D=3, **kwargs)
    reconstruct_angles(binning=5, D=4, **kwargs)
    reconstruct_angles(binning=5, D=5, **kwargs)
    kwargs['randomize_binning'] = True
    reconstruct_angles(binning=1, D=1, **kwargs)
    reconstruct_angles(binning=1, D=2, **kwargs)
    reconstruct_angles(binning=1, D=3, **kwargs)
    reconstruct_angles(binning=1, D=4, **kwargs)
    reconstruct_angles(binning=1, D=5, **kwargs)
    reconstruct_angles(binning=2.5, D=1, **kwargs)
    reconstruct_angles(binning=2.5, D=2, **kwargs)
    reconstruct_angles(binning=2.5, D=3, **kwargs)
    reconstruct_angles(binning=2.5, D=4, **kwargs)
    reconstruct_angles(binning=2.5, D=5, **kwargs)
    reconstruct_angles(binning=5, D=1, **kwargs)
    reconstruct_angles(binning=5, D=2, **kwargs)
    reconstruct_angles(binning=5, D=3, **kwargs)
    reconstruct_angles(binning=5, D=4, **kwargs)
    reconstruct_angles(binning=5, D=5, **kwargs)

def reconstruct_angles(data, tablename, dest, THETA, D, binning=False,
                       randomize_binning=False, N=None):
    """Reconstruct angles from simulation for minimum particle density"""

    match = re.search('_size([0-9]+)', tablename)
    if match:
        R = int(match.group(1))
    else:
        R = 10

    table = data.getNode('/analysis', tablename)
    dst_row = dest.row
    for event in table[:N]:
        if min(event['n1'], event['n3'], event['n4']) >= D:
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
            alpha = event['alpha']

            if not isnan(theta) and not isnan(phi):
                ang_dist = arccos(sin(theta) * sin(THETA) *
                                  cos(phi - -alpha) + cos(theta) *
                                  cos(THETA))

                dst_row['r'] = event['r']
                dst_row['phi'] = event['phi']
                dst_row['alpha'] = alpha
                dst_row['t1'] = event['t1']
                dst_row['t2'] = event['t2']
                dst_row['t3'] = event['t3']
                dst_row['t4'] = event['t4']
                dst_row['n1'] = event['n1']
                dst_row['n2'] = event['n2']
                dst_row['n3'] = event['n3']
                dst_row['n4'] = event['n4']
                dst_row['sim_theta'] = THETA
                dst_row['sim_phi'] = -alpha
                dst_row['r_theta'] = theta
                dst_row['r_phi'] = phi
                dst_row['D'] = min(event['n1'], event['n3'], event['n4'])
                dst_row['size'] = R
                if binning is False:
                    bin_size = 0
                else:
                    bin_size = binning
                dst_row['bin'] = bin_size
                dst_row['bin_r'] = randomize_binning
                dst_row.append()
    dest.flush()

def do_reconstruction_plots(data, tablename):
    """Make plots based upon earlier reconstructions"""

    table = data.getNode('/reconstructions', tablename)

    figure()
    x, y, y2 = [], [], []
    for D in range(1, 6):
        x.append(D)
        events = table.readWhere(
            '(D==%d) & (sim_theta==%.40f) & (size==10) & (bin==0)' % 
            (D, float32(pi / 8)))
        print len(events),
        errors = events['sim_theta'] - events['r_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['sim_phi'] - events['r_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(x, rad2deg(y), '.-', label="Theta")
    plot(x, rad2deg(y2), '.-', label="Phi")
    xlabel("Minimum number of particles")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title("Uncertainty as a function of number of particles")
    figtext(.65, .8, "Azimuthal angle: 22.5 degrees",
            horizontalalignment='right')
    legend()
    savefig('plots/auto-results-MIP.pdf')
    print

    figure()
    x, y, y2 = [], [], []
    for THETA in [0, deg2rad(5), pi / 8, deg2rad(35)]:
        x.append(THETA)
        events = table.readWhere(
            '(D>=2) & (sim_theta==%.40f) & (size==10) & (bin==0)' % 
            float32(THETA))
        print len(events),
        errors = events['sim_theta'] - events['r_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['sim_phi'] - events['r_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(rad2deg(x), rad2deg(y), '.-', label="Theta")
    # Azimuthal angle undefined for zenith = 0
    plot(rad2deg(x[1:]), rad2deg(y2[1:]), '.-', label="Phi")
    xlabel("Shower zenith angle (degrees)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title("Uncertainty as a function of shower zenith angle")
    figtext(.65, .8, "Number of particles at least 2",
            horizontalalignment='right')
    legend()
    ylim(ymin=0)
    savefig('plots/auto-results-zenith.pdf')
    print

    figure()
    x, y, y2 = [], [], []
    for size in [5, 10, 20]:
        x.append(size)
        events = table.readWhere(
            '(D>=2) & (sim_theta==%.40f) & (size==%d) & (bin==0)' %
            (float32(pi/ 8), size))
        print len(events),
        errors = events['sim_theta'] - events['r_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['sim_phi'] - events['r_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(x, rad2deg(y), '.-', label="Theta")
    plot(x, rad2deg(y2), '.-', label="Phi")
    xlabel("Station size (m)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title("Uncertainty as a function of station size")
    figtext(.65, .8, "Number of particles at least 2\n"
            "Azimuthal angle: 22.5 degrees", horizontalalignment='right')
    legend()
    savefig('plots/auto-results-size.pdf')
    print

    figure()
    x, y, y2 = [], [], []
    for bin_size in [0, 1, 2.5, 5]:
        x.append(bin_size)
        events = table.readWhere(
            '(D>=2) & (sim_theta==%.40f) & (size==10) & (bin==%.40f) & '
            '(bin_r==False)' %
            (float32(pi/ 8), bin_size))
        print len(events),
        errors = events['sim_theta'] - events['r_theta']
        # Make sure -pi < errors < pi
        errors = (errors + pi) % (2 * pi) - pi
        errors2 = events['sim_phi'] - events['r_phi']
        # Make sure -pi < errors2 < pi
        errors2 = (errors2 + pi) % (2 * pi) - pi
        y.append(std(errors))
        y2.append(std(errors2))
    plot(x, rad2deg(y), '.-', label="Theta")
    plot(x, rad2deg(y2), '.-', label="Phi")
    xlabel("Bin size (ns)")
    ylabel("Uncertainty in angle reconstruction (deg)")
    title("Uncertainty as a function of bin size")
    figtext(.65, .8, "Number of particles at least 2\n"
            "Azimuthal angle: 22.5 degrees", horizontalalignment='right')
    legend(loc='best')
    ylim(ymin=0)
    savefig('plots/auto-results-binsize.pdf')
    print


if __name__ == '__main__':
    # invalid values in arcsin will be ignored (nan handles the situation
    # quite well)
    np.seterr(invalid='ignore')

    try:
        data
    except NameError:
        data = tables.openFile('data-e15.h5', 'a')

    #do_full_reconstruction(data, 'full')
    do_reconstruction_plots(data, 'full')
