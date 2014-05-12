import tables

from pylab import *
from numpy import *

from artist import GraphArtist


def get_front_arrival_time(sim, R, dR, theta):
    query = '(R - dR <= core_distance) & (core_distance < R + dR)'
    c = 3e-1 # m / ns

    t_list = []
    for shower in sim:
        particles = shower.leptons.read_where(query)
        x = particles[:]['x']
        t = particles[:]['arrival_time']

        dt = x * sin(theta) / c
        t += dt

        t_list.extend(t)

    return array(t_list)


def monte_carlo_timings(n, bins, size):
    x0 = min(bins)
    x1 = max(bins)
    y0 = 0
    y1 = max(n)

    t_list = []
    while len(t_list) < size:
        x = random.uniform(x0, x1)
        y = random.uniform(y0, y1)
        idx = bins.searchsorted(x) - 1
        if y <= n[idx]:
            t_list.append(x)

    return t_list


@vectorize
def my_std_t(data, N):
    sim = data.root.showers.E_1PeV.zenith_22_5
    t = get_front_arrival_time(sim, 30, 5, pi / 8)
    n, bins = histogram(t, bins=linspace(0, 50, 401))
    mct = monte_carlo_timings(n, bins, 10000)
    print "Monte Carlo:", N

    mint_list = []
    i = 0
    while i < len(mct):
        try:
            values = mct[i:i + N]
        except IndexError:
            break
        if len(values) == N:
            mint_list.append(min(values))
        i += N
    return median(mint_list)


def my_std_t_for_R(data, N_list, R_list):
    sim = data.root.showers.E_1PeV.zenith_22_5

    value_list = []
    for N, R in zip(N_list, R_list):
        t = get_front_arrival_time(sim, R, 5, pi / 8)
        n, bins = histogram(t, bins=linspace(0, 50, 401))
        mct = monte_carlo_timings(n, bins, 10000)
        print "Monte Carlo:", N

        mint_list = []
        i = 0
        while i < len(mct):
            try:
                values = mct[i:i + N]
            except IndexError:
                break
            if len(values) == N:
                mint_list.append(min(values))
            i += N
        value_list.append(median(mint_list))
    return array(value_list)


def my_t_draw_something(data, N, num_events):
    sim = data.root.showers.E_1PeV.zenith_22_5
    t = get_front_arrival_time(sim, 20, 5, pi / 8)
    n, bins = histogram(t, bins=linspace(0, 50, 201))
    mct = monte_carlo_timings(n, bins, num_events * N)
    print "Monte Carlo:", N

    mint_list = []
    i = 0
    while i < len(mct):
        try:
            values = mct[i:i + N]
        except IndexError:
            break
        if len(values) == N:
            mint_list.append(min(values))
        i += N
    return mint_list


def plot_R():
    graph = GraphArtist(width=r'.45\linewidth')

    n, bins, patches = hist(data.root.simulations.E_1PeV.zenith_22_5.shower_0.coincidences.col('r'), bins=100, histtype='step')
    graph.histogram(n, bins, linestyle='black!50')

    shower = data.root.simulations.E_1PeV.zenith_22_5.shower_0
    ids = shower.observables.get_where_list('(n1 >= 1) & (n3 >= 1) & (n4 >= 1)')
    R = shower.coincidences.read_coordinates(ids, field='r')
    n, bins, patches = hist(R, bins=100, histtype='step')
    graph.histogram(n, bins)

    xlabel("Core distance [m]")
    ylabel("Number of events")

    print "mean", mean(R)
    print "median", median(R)

    graph.set_xlabel(r"Core distance [\si{\meter}]")
    graph.set_ylabel("Number of events")
    graph.set_xlimits(min=0)
    graph.set_ylimits(min=0)
    graph.save('plots/SIM-R')


def plot_arrival_times():
    graph = GraphArtist()

    figure()
    sim = data.root.showers.E_1PeV.zenith_22_5
    t = get_front_arrival_time(sim, 20, 5, pi / 8)
    n, bins = histogram(t, bins=linspace(0, 50, 201))
    mct = monte_carlo_timings(n, bins, 100000)
    n, bins, patches = hist(mct, bins=linspace(0, 20, 101), histtype='step')
    graph.histogram(n, bins, linestyle='black!50')

    mint = my_t_draw_something(data, 2, 100000)
    n, bins, patches = hist(mint, bins=linspace(0, 20, 101), histtype='step')
    graph.histogram(n, bins)

    xlabel("Arrival time [ns]")
    ylabel("Number of events")

    graph.set_xlabel(r"Arrival time [\si{\nano\second}]")
    graph.set_ylabel("Number of events")
    graph.set_xlimits(0, 20)
    graph.set_ylimits(min=0)
    graph.save('plots/SIM-T')

    print median(t), median(mct), median(mint)


if __name__ == '__main__':
    if not 'data' in globals():
        data = tables.open_file('master-ch4v2.h5')

    plot_R()
    plot_arrival_times()
