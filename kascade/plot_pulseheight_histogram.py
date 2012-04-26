import tables

from pylab import *

from artist import GraphArtist


def plot_pulseheight_histogram(data):
    events = data.root.hisparc.cluster_kascade.station_601.events
    ph = events.col('pulseheights')

    clf()
    n, bins, patches = hist(ph[:, 0], bins=arange(0, 2001, 10),
                            histtype='step')
    yscale('log')
    ylim(ymin=1e1)

    graph = GraphArtist('semilogy', width=r'.5\linewidth')
    graph.histogram(n, bins)
    graph.set_xlabel(r"Pulseheight [\adc{}]")
    graph.set_ylabel(r"Number of events")
    graph.save("plots/plot_pulseheight_histogram")


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.openFile('kascade.h5')

    plot_pulseheight_histogram(data)
