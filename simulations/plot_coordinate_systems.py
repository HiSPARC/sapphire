from pylab import *

import utils


def plot_coordinate_systems():
    figure()
    suptitle('positions and orientations')

    x, y, alpha = generate_random_coordinates_in_circle(10, 1000)
    xp, yp, alphap = transform_coordinates(x, y, alpha)

    subplot('121', aspect='equal')
    draw_coordinates(x, y, alpha)
    title('shower-centered coordinates')

    subplot('122', aspect='equal')
    draw_coordinates(xp, yp, alphap)
    title('cluster-centered coordinates')

    utils.saveplot()


def generate_coordinates_in_rectangle():
    x, y, alpha = [], [], []
    for i in linspace(-10, 10, 5):
        for j in linspace(-10, 10, 5):
            for k in linspace(-pi, pi, 11):
                x.append(i)
                y.append(j)
                alpha.append(k)

    return array(x), array(y), array(alpha)


def transform_coordinates(x, y, alpha):
    xp, yp, alphap = [], [], []

    for u, v, w in zip(x, y, alpha):
        r = sqrt(u ** 2 + v ** 2)
        phi = arctan2(v, u)

        phi += pi - w

        # shower azimuth angle, if it was initially zero
        alpha = -w

        xp.append(r * cos(phi))
        yp.append(r * sin(phi))
        alphap.append(alpha)

    return array(xp), array(yp), array(alphap)


def draw_coordinates(x, y, alpha):
    for u, v, w in zip(x, y, alpha):
        scatter(u, v, c='black', s=4)
        s = sin(w)
        t = cos(w)
        plot([u, u + s], [v, v + t], c='black')


def plot_coordinate_density():
    figure()
    suptitle('densities')

    x, y, alpha = generate_random_coordinates_in_circle(10, 100000)
    xp, yp, alphap = transform_coordinates(x, y, alpha)

    subplot('121', aspect='equal')
    draw_coordinate_density(x, y)
    title('shower-centered coordinates')

    subplot('122', aspect='equal')
    draw_coordinate_density(xp, yp)
    title('cluster-centered coordinates')

    utils.saveplot()


def generate_random_coordinates_in_rectangle(N=100):
    x, y = uniform(-10, 10, (2, N))
    alpha = uniform(-pi, pi, N)

    return x, y, alpha


def generate_random_coordinates_in_circle(R, N=100):
    x, y, alpha = [], [], []

    while len(x) < N:
        u, v = uniform(-10, 10, 2)

        if u ** 2 + v ** 2 <= R ** 2:
            x.append(u)
            y.append(v)
            alpha.append(uniform(-pi, pi))

    return x, y, alpha


def draw_coordinate_density(x, y):
    H, xedges, yedges = histogram2d(x, y, bins=100)
    contourf(H.T, 10, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])


if __name__ == '__main__':
    plot_coordinate_systems()
    plot_coordinate_density()
