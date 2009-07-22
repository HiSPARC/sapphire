""" Clean up coincidences

    The following snippet cleans up the HiSPARC / KASCADE coincidences by
    removing events which are far apart.  That is, events which are only
    considered nearest neighbours because the HiSPARC detector was down for
    some time.

"""
if __name__ == '__main__':
    coincidences = data.root.coincidences.events
    l = []
    for i in range(len(coincidences)):
        c = coincidences[i]
        t0 = c['hisparc_ext_timestamp']
        t1 = c['kascade_ext_timestamp']
        dt = abs(t0 - t1)

        if dt > int(14e9):
            l.append(i)
