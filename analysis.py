import tables
import progressbar as pb
import inspect

import sys
sys.path.append('src/')

from clusters import SimpleCluster


def whosparent():
    """Return parent function name of caller"""

    return inspect.stack()[2][3]

def saveplot(suffix=''):
    """Save a plot using caller's name"""

    if suffix:
        suffix = '-' + suffix
    savefig('plots/%s%s%s.pdf' % (whosparent(), suffix, __suffix))

def title(text):
    pylab.title(text + '\n(%s)' % __suffix)

mylog = vectorize(lambda x: log10(x) if x > 0 else 0)


def main(group, suffix):
    global __suffix
    __suffix = suffix

    hist_Rcluster(group)
    hist_Rstation(group)
    hist_mean_Nparticles(group)
    hist2d_stations(group)
    hist2d_clusters(group)
    scatter_cores_altcs(group)

def hist_Rcluster(group):
    figure()
    coincidences = group.coincidences

    for Ntrig in range(1, 5):
        sel = coincidences.readWhere('N == Ntrig')
        hist(sel['r'], bins=200, histtype='step', label="N=%d" % Ntrig)

    xlabel("R_cluster [m]")
    ylabel("Count")
    title("Core distance to cluster center")
    legend(loc='best')
    saveplot()

def hist_Rstation(group):
    figure()
    coincidences = group.coincidences
    c_index = group.c_index
    observables = group.observables

    for Ntrig in range(1, 5):
        R = []
        sel = coincidences.readWhere('N == Ntrig')
        progress = pb.ProgressBar(widgets=[pb.Percentage(), pb.Bar(),
                                           pb.ETA()])
        for event_id in progress(sel['id']):
            R.extend([observables[u]['r'] for u in c_index[event_id]])
        hist(R, bins=200, histtype='step', label="N=%d" % Ntrig)

    xlabel("R_station [m]")
    ylabel("Count")
    title("Core distance to station center")
    legend(loc='best')
    saveplot()

def hist_mean_Nparticles(group):
    figure()
    coincidences = group.coincidences
    c_index = group.c_index
    observables = group.observables

    sel = coincidences.readWhere('N >= 3')
    mean_n = []
    for event_id in sel['id']:
        mean_n.extend([mean([observables[u][v] for v in ['n1', 'n2', 'n3',
                                                         'n4']]) for
                       u in c_index[event_id]])
    hist(mean_n, bins=arange(0, 10, .25), histtype='step')

    xlabel("Mean N_particles")
    ylabel("Count")
    title("Mean N_particles per detector (N_station >= 3)")
    saveplot()

def hist2d_stations(group):
    figure()
    observables = group.observables

    sel = observables.readWhere('N >= 2')
    H, xedges, yedges = histogram2d(sel['x'], sel['y'], bins=100)
    contourf(H.T, 20, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    xlabel("Distance [m]")
    ylabel("Distance [m]")
    title("Occurences of triggered stations (N >= 2)")
    saveplot()

def hist2d_clusters(group):
    figure()
    coincidences = group.coincidences

    sel = coincidences.readWhere('N >= 3')
    H, xedges, yedges = histogram2d(sel['x'], sel['y'], bins=100)
    contourf(H.T, 20, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
    xlabel("Distance [m]")
    ylabel("Distance [m]")
    title("Occurences of triggered clusters (N >= 3)")
    saveplot()

def scatter_cores_altcs(group):
    figure()
    coincidences = group.coincidences

    sel = coincidences.readWhere('N >= 3')
    xp = []
    yp = []
    for event in sel:
        rp = event['r']
        phip = event['phi'] + pi
        phip -= event['alpha']
        phip = (phip + pi) % (2 * pi) - pi
        xp.append(rp * cos(phip))
        yp.append(rp * sin(phip))
    plot(xp, yp, ',')

    draw_simplecluster(0, 0, 0, 'r', draw_pos=True)

    xlabel("Distance [m]")
    ylabel("Distance [m]")
    title("Shower core positions (N >= 3)")
    saveplot()

def draw_simplecluster(r, phi, alpha, spec='r', draw_pos=False):
    cluster = SimpleCluster(size=150)
    for station in cluster.stations:
        x, y, beta = station.get_coordinates(r, phi, alpha)
        for detector in station.detectors:
            c = detector.get_corners(x, y, beta)
            cx, cy = zip(*c)
            fill(cx, cy, spec, ec='none')
        if draw_pos:
            scatter(x, y, c=spec)


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('data-e15-S150.h5', 'r')

    sim = data.root.simulations.E_1PeV.zenith_0
    main(sim, 'E-1PeV-S150m')
