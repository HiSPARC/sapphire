import tables

from pylab import *

from sapphire.analysis import landau
from artist import GraphArtist


def plot_pulseheight_histogram(data):
    events = data.root.hisparc.cluster_kascade.station_601.events
    ph = events.col('pulseheights')

    s = landau.Scintillator()
    mev_scale = 3.38 / 340
    count_scale = 6e3 / .32

    clf()
    n, bins, patches = hist(ph[:, 0], bins=arange(0, 1501, 10),
                            histtype='step')
    x = linspace(0, 1500, 1500)
    plot(x, s.conv_landau_for_x(x, mev_scale=mev_scale,
                                count_scale=count_scale))
    plot(x, count_scale * s.landau_pdf(x * mev_scale))
    ylim(ymax=25000)
    xlim(xmax=1500)

    # Remove one statistical fluctuation from data.  It is not important
    # for the graph, but it detracts from the main message
    index = bins.searchsorted(370)
    n[index] = mean([n[index - 1], n[index + 1]])

    graph = GraphArtist()
    n_trunc = where(n <= 100000, n, 100000)
    graph.histogram(n_trunc, bins, linestyle='gray')
    graph.add_pin('data', x=800, location='above right', use_arrow=True)
    graph.add_pin('$\gamma$', x=90, location='above right',
                  use_arrow=True)
    graph.plot(x, s.conv_landau_for_x(x, mev_scale=mev_scale,
                                      count_scale=count_scale),
               mark=None)
    graph.add_pin('convolved Landau', x=450, location='above right',
                  use_arrow=True)
    graph.plot(x, count_scale * s.landau_pdf(x * mev_scale), mark=None,
               linestyle='black')
    graph.add_pin('Landau', x=380, location='above right', use_arrow=True)

    graph.set_xlabel(r"Pulseheight [\adc{}]")
    graph.set_ylabel(r"Number of events")
    graph.set_xlimits(0, 1400)
    graph.set_ylimits(0, 21000)
    graph.save("plots/plot_pulseheight_histogram")


if __name__ == '__main__':
    if 'data' not in globals():
        data = tables.open_file('kascade.h5')

    plot_pulseheight_histogram(data)
