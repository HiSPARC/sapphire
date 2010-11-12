import zlib
import tables
from numpy import *
from pylab import *

from IPython.Shell import IPShellEmbed
ipshell = IPShellEmbed()

ADC_THRESHOLD = 20
ADC_TIME_PER_SAMPLE = 2.5e-9

DETECTORS = [(65., 15.05, 'UD'), (65., 20.82, 'UD'), (70., 23.71, 'LR'),
             (60., 23.71, 'LR')]


def process_traces(events, traces_table, limit=None):
    """Process traces to yield pulse timing information"""

    result = []
    result2 = []
    for event in events[:limit]:
        trace = get_traces(traces_table, event)
        timings, timings2 = zip(*[reconstruct_time_from_trace(x) for x in
                                  trace])
        result.append(timings)
        result2.append(timings2)
    return 1e9 * array(result), 1e9 * array(result2)

def get_traces(traces_table, event):
    """Retrieve traces from table and reconstruct them"""

    if type(event) != list:
        idxs = event['traces']
    else:
        idxs = event

    traces = []
    for idx in idxs:
        trace = zlib.decompress(traces_table[idx]).split(',')
        if trace[-1] == '':
            del trace[-1]
        trace = array([int(x) for x in trace])
        traces.append(trace)
    return traces

def reconstruct_time_from_trace(trace):
    """Reconstruct time of measurement from a trace"""

    t = trace[:100]
    baseline = mean(t)
    stdev = std(t)

    trace = trace - baseline
    threshold = ADC_THRESHOLD

    value = nan
    for i, t in enumerate(trace):
        if t >= threshold:
            value = i
            break

    # Better value, interpolation
    if not isnan(value):
        x0, x1 = i - 1, i
        y0, y1 = trace[x0], trace[x1]
        v2 = 1. * (threshold - y0) / (y1 - y0) + x0
    else:
        v2 = nan

    return value * ADC_TIME_PER_SAMPLE, v2 * ADC_TIME_PER_SAMPLE

def plot_traces(data, station, event_id):
    """Plot traces from a single event"""

    parent = data.getNode('/', station)
    events = data.getNode(parent, 'events')
    traces_table = data.getNode(parent, 'blobs')

    figure()
    traces = get_traces(traces_table, events[event_id - 1])
    x = arange(0, len(traces[0])) * 2.5e-3
    for n, trace in enumerate(traces, 1):
        plot(x, trace * -.57 + 114, label="Scint %d" % n)

    xlabel("Time (us)")
    ylabel("Voltage (mV)")
    title("Traces from event %d (station %s)" % (event_id, station))
    legend(loc='best')

def frac_bins(low, high, binsize, nbins=1):
    binsize = binsize * nbins
    low = low - .5 * binsize
    high = ceil((high - low) / binsize) * binsize + low + .5 * binsize
    return arange(low, high, binsize)

def get_timing_data(data):
    """Get timing data from data file or analyze it now"""

    try:
        timing_data = data.root.analysis.timing_data.read()
    except tables.NoSuchNodeError:
        print 'no timing data analyzed, analyzing now...'
        events = data.getNode('/kascade', 'coincidences')
        traces_table = data.getNode('/',
                            'hisparc/cluster_kascade/station_601/blobs')
        timing_data = process_traces(events, traces_table)
        print "storing results for future use"
        if not data.__contains__('/analysis'):
            data.createGroup('/', 'analysis', "Analyzed data")
        data.createArray('/analysis', 'timing_data', timing_data,
                         "Start-of-pulse timings derived from traces")

    return timing_data

def plot_reconstructed_angles(events, timing_data, D=2., s='', shifts=None,
                              randomize=False):
    theta_list = []
    phi_list = []

    k_theta_list = []
    k_phi_list = []

    NT, NS, NF = 0, 0, 0
    for event, timings in zip(events[:], timing_data):
        if shifts:
            timings += shifts

        ph1, ph2, ph3, ph4 = event['pulseheights']
        if min([ph1 / 380., ph3 / 410., ph4 / 380.]) >= D:
            NT += 1
            if randomize is True:
                timings = [x + uniform(0, 2.5) for x in timings]
            theta, phi = reconstruct_angle(timings)
            if not isnan(theta) and not isnan(phi):
                NS += 1
                k_theta = event['k_zenith']
                k_phi = event['k_azimuth']
                theta_list.append(theta)
                phi_list.append(phi)
                k_theta_list.append(k_theta)
                k_phi_list.append(k_phi)
            else:
                NF += 1

    thetas = array(theta_list)
    k_thetas = array(k_theta_list)
    phis = array(phi_list)
    k_phis = array(k_phi_list)
    k_phis_orig = k_phis
    k_phis = -(k_phis + deg2rad(75)) % (2 * pi) - pi

    angular_dists = arccos(sin(thetas) * sin(k_thetas) *
                           cos(phis - k_phis) +
                           cos(thetas) * cos(k_thetas))


    figure()
    plot(rad2deg(k_thetas), rad2deg(thetas), '.', ms=1.)
    title("Theta angle reconstruction (min. D >= %.1f)" % D)
    xlabel("KASCADE theta angle")
    xlim(0, 40)
    ylabel("HiSPARC theta angle")
    ylim(0, 90)
    savefig("plots/auto-theta-D%d%s.pdf" % (D, s))

    figure()
    plot(rad2deg(k_phis), rad2deg(phis), '.', ms=1.)
    title("Phi angle reconstruction (min. D >= %.1f)" % D)
    xlabel("KASCADE phi angle")
    xlim(-180, 180)
    ylabel("HiSPARC phi angle")
    ylim(-180, 180)
    savefig("plots/auto-phi-D%d%s.pdf" % (D, s))

    figure()
    plot(rad2deg(k_phis_orig), rad2deg(phis), '.', ms=1.)
    title("Phi angle reconstruction (min. D >= %.1f)" % D)
    xlabel("KASCADE phi angle")
    xlim(0, 360)
    ylabel("HiSPARC phi angle")
    ylim(-180, 180)
    savefig("plots/auto-phiorig-D%d%s.pdf" % (D, s))

    figure()
    hist(rad2deg(thetas - k_thetas), bins=200, histtype='step')
    title("Theta angle reconstruction accuracy (min. D >= %.1f)" % D)
    xlabel("HiSPARC - KASCADE theta angle")
    ylabel("count")
    savefig("plots/auto-thetahist-D%d%s.pdf" % (D, s))

    figure()
    hist(rad2deg(phis - k_phis), bins=200, histtype='step')
    title("Phi angle reconstruction accuracy (min. D >= %.1f)" % D)
    xlabel("HiSPARC - KASCADE phi angle")
    ylabel("count")
    savefig("plots/auto-phihist-D%d%s.pdf" % (D, s))

    figure()
    n, bins, patches = hist(rad2deg(angular_dists),
                            bins=linspace(0, 120, 200), histtype='step')
    res = get_resolution(n, bins)
    #axvspan(xmin=0, xmax=res, color='blue', alpha=.2)
    axvline(x=res, color='blue', alpha=.2)
    title("Angular distance of HiSPARC and KASCADE angles (D >= %d)" % D)
    xlabel("Angular distance (degrees)")
    ylabel("Count")
    figtext(.65, .8, "resolution: %.2f deg\nfailed: %5.1f %%" %
            (res, 100. * NF / NT))
    savefig("plots/auto-angdist-D%d%s.pdf" % (D, s))

    figure()
    hist(rad2deg(thetas), bins=200, histtype='step')
    title("KASCADE zenith angles")
    xlabel("zenith angle (degrees)")
    ylabel("count")

    print
    print
    print "Angle reconstruction (D >= %d)" % D
    print "Shifts: %r" % shifts
    print "Tag: %s" % s
    print "Total of %d events" % len(events)
    print "Total number reconstructions:        %6d" % NT
    print "Number of succesful reconstructions: %6d (%5.1f %%)" % \
        (NS, 100. * NS / NT)
    print "Number of failed reconstructions:    %6d (%5.1f %%)" % \
        (NF, 100. * NF / NT)
    print "Angle resolution (68%% integrated): %.2f degrees" % res
    print

def reconstruct_angle(event_timing):
    """Reconstruct angles from a single event"""

    c = 3.00e+8

    #dt1 = event['t1'] - event['t3']
    #dt2 = event['t1'] - event['t4']

    dt1 = event_timing[0] - event_timing[2]
    dt2 = event_timing[0] - event_timing[3]

    r1, phi1 = calc_r_phi(1, 3)
    r2, phi2 = calc_r_phi(1, 4)

    phi = arctan2((dt2 * r1 * cos(phi1) - dt1 * r2 * cos(phi2)),
                  (dt2 * r1 * sin(phi1) - dt1 * r2 * sin(phi2)) * -1)
    theta = arcsin(c * dt1 * 1e-9 / (r1 * cos(phi - phi1)))

    return theta, phi

def calc_r_phi(s1, s2):
    """Calculate (r, phi) for two detectors s1 and s2"""

    x1, y1 = DETECTORS[s1 - 1][:2]
    x2, y2 = DETECTORS[s2 - 1][:2]

    return sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2), \
           arctan2((y2 - y1), (x2 - x1))

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


if __name__ == '__main__':
    # invalid values in arcsin will be ignored (nan handles the situation
    # quite well)
    np.seterr(invalid='ignore')

    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'a')

    try:
        timing_data
        timing_data_linear
    except NameError:
        timing_data = data.root.analysis_new.timing_data.read()
        timing_data_linear = data.root.analysis_new.timing_data_linear.read()

    events = data.root.kascade_new.coincidences
    print "Reconstructing angles..."
    plot_reconstructed_angles(events, timing_data, 1., randomize=False)
    plot_reconstructed_angles(events, timing_data, 2., randomize=False)
    plot_reconstructed_angles(events, timing_data, 4., randomize=False)

    plot_reconstructed_angles(events, timing_data_linear, 1., '-linear',
                              randomize=False)
    plot_reconstructed_angles(events, timing_data_linear, 2., '-linear',
                              randomize=False)
    plot_reconstructed_angles(events, timing_data_linear, 4., '-linear',
                              randomize=False)

    plot_reconstructed_angles(events, timing_data, 1., '-shifts',
                              shifts=[.25, 0., 1.17, -.21],
                              randomize=False)
    plot_reconstructed_angles(events, timing_data, 2., '-shifts',
                              shifts=[.25, 0., 1.17, -.21],
                              randomize=False)
    plot_reconstructed_angles(events, timing_data, 4., '-shifts',
                              shifts=[.25, 0., 1.17, -.21],
                              randomize=False)

    plot_reconstructed_angles(events, timing_data_linear, 1.,
                              '-linear-shifts',
                              shifts=[.26, 0., 1.20, -.19],
                              randomize=False)
    plot_reconstructed_angles(events, timing_data_linear, 2.,
                              '-linear-shifts',
                              shifts=[.26, 0., 1.20, -.19],
                              randomize=False)
    plot_reconstructed_angles(events, timing_data_linear, 4.,
                              '-linear-shifts',
                              shifts=[.26, 0., 1.20, -.19],
                              randomize=False)
