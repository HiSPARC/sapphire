import tables

from pylab import *

from hisparc.analysis.traces import get_traces
from artist import GraphArtist


def plot_trace(station_group, idx):
    events = station_group.events
    blobs = station_group.blobs

    traces_idx = events[idx]['traces']
    traces = get_traces(blobs, traces_idx)
    traces = array(traces)
    x = arange(traces.shape[1])
    x *= 2.5

    clf()
    plot(x, traces.T)
    xlim(0, 200)

    #line_styles = ['solid', 'dashed', 'dotted', 'dashdotted']
    line_styles = ['black', 'black!80', 'black!60', 'black!40']
    styles = (u for u in line_styles)

    graph = GraphArtist(width=r'.5\linewidth')
    for trace in traces:
        graph.plot(x, trace / 1000, mark=None, linestyle=styles.next())
    graph.set_xlabel(r"Time [\si{\nano\second}]")
    graph.set_ylabel(r"Signal [\si{\volt}]")
    graph.set_xlimits(0, 200)
    graph.save('plots/traces')


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.open_file('master.h5')

    station = data.root.s501
    ph = station.events.col('pulseheights')
    phmax = ph.max(1)
    idxs = (phmax > 2000).nonzero()[0]

    plot_trace(station, idxs[11])
