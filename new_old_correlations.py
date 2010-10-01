import tables

from numpy import *
from pylab import *


def mylog(x):
    if x >= 10.:
        return log10(x)
    elif 0 < x < 10:
        return x / 10.
    else:
        return 0.
vlog = vectorize(mylog)

def plot_new_old_correlations(data):
    c = data.root.kascade.coincidences
    cn = data.root.kascade_new.coincidences

    figure()
    H, xedges, yedges = histogram2d(c[:]['k_dens_e'][:,0],
                                    c[:]['pulseheights'][:,0],
                                    range=[[0, 10], [150, 3000]], bins=100)
    contourf([mean([u, v]) for u, v in zip(xedges[:-1], xedges[1:])],
             [mean([u, v]) for u, v in zip(yedges[:-1], yedges[1:])],
             vlog(H.T), 20)
    colorbar()
    title("HiSPARC / KASCADE correlation (detector 1) (OLD SET)")
    xlabel("KASCADE electron density (m^{-2})")
    ylabel("HiSPARC pulseheight (ADC)")

    figure()
    H, xedges, yedges = histogram2d(cn[:]['k_dens_e'][:,0],
                                    cn[:]['pulseheights'][:,0],
                                    range=[[0, 10], [150, 3000]], bins=100)
    contourf([mean([u, v]) for u, v in zip(xedges[:-1], xedges[1:])],
             [mean([u, v]) for u, v in zip(yedges[:-1], yedges[1:])],
             vlog(H.T), 20)
    colorbar()
    title("HiSPARC / KASCADE correlation (detector 1) (NEW SET)")
    xlabel("KASCADE electron density (m^{-2})")
    ylabel("HiSPARC pulseheight (ADC)")


if __name__ == '__main__':
    try:
        data
    except NameError:
        data = tables.openFile('kascade.h5', 'a')

    plot_new_old_correlations(data)
