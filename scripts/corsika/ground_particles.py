#!/usr/bin/env python

import matplotlib.pyplot as plt

from sapphire import corsika


def plot_ground(x, y, eventheader, title='Ground particles'):
    size = 200.
    plt.figure(figsize=(9, 9))
    plt.scatter(x, y, c='r', s=2., edgecolors='none')
    plt.axis('equal')
    plt.axis([-size, size, -size, size])
    plt.xlabel('x [m]')
    plt.ylabel('y [m]')
    plt.title(title)
#     plt.savefig('/Users/arne/Dropbox/hisparc/Plots/Corsika/scatter_em.pdf')
    plt.show()
    return


def main():
    corsika_data = corsika.CorsikaFile('../../sapphire/tests/corsika/DAT000000')
    corsika_data.check()

    for event in corsika_data.get_events():
        x = []
        y = []
        for particle in event.get_particles():
            if particle.is_detectable:
                x.append(particle.x)
                y.append(particle.y)
        event_header = event.get_header()
        title = ('Ground particles, Primary: %s, E = %1.0e eV' %
                 (event_header.particle, event_header.energy))
        plot_ground(x, y, event_header, title=title)


if __name__ == '__main__':
    main()
