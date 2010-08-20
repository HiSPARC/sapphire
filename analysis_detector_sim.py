import tables
from itertools import combinations

from pylab import *


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


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('simulation-e15.h5', 'r')

    plot_all_ring_timings(data)
