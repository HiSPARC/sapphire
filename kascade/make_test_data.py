"""Simple utility to make a truncated copy of the full KASCADE data"""

import os

import tables


N = 100000


if __name__ == '__main__':
    data = tables.openFile('kascade.h5', 'r')
    data2 = tables.openFile('kascade-small.h5', 'w')

    group = '/hisparc/cluster_kascade/station_601'
    data2.createGroup(*os.path.split(group), createparents=True)

    events = data.getNode(group, 'events')
    events2 = data2.createTable(group, 'events', events.description)
    blobs = data.getNode(group, 'blobs')
    blobs2 = data2.createVLArray(group, 'blobs', blobs.atom)

    for event in events.iterrows(stop=N):
        events2.append([event[:]])
        for idx in event['traces']:
            if idx != -1:
                blobs2.append(blobs[idx])
