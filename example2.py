#!/usr/bin/env python

import matplotlib.pyplot as plt
import corsika


def plot_ground(x, y, eventheader, title='Ground particles'):
    size = 200.
    plt.figure(figsize=(9, 9))
    plt.scatter(x, y, c='r', s=2., edgecolors='none')
    plt.axis('equal')
    plt.axis([-size, size, -size, size])
    plt.xlabel('x in meters')
    plt.ylabel('y in meters')
    energy = r'E = %s eV' % eventheader.fEnergy

    plt.title(title)
    plt.savefig('/Users/arne/Dropbox/hisparc/Plots/Corsika/scatter_em.pdf')
    plt.show()
    return


def main():
    corsika_data = corsika.CorsikaFile('/Users/arne/Datastore/CORSIKA/DAT000001')
    corsika_data.Check()

    for event in corsika_data.GetEvents():
        x = []
        y = []
        for particle in event.GetParticles():
            if particle.IsEM():
                x.append(particle.fX)
                y.append(particle.fY)
        EventHeader = event.GetHeader()
        title = 'Ground EM particles, %s shower' % EventHeader.fEnergy
        plot_ground(x, y, EventHeader, title=title)


if __name__ == '__main__':
    main()
