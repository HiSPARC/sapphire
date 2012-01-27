from numpy import *
from pylab import *

def contour_histogram2d(x, y, bins=50, levels=10, log=False):
    mylog = vectorize(lambda x: log10(x) if x > 0 else 0.)

    H, xedges, yedges = histogram2d(x, y, bins)
    if log:
        H = mylog(H)
    x = (xedges[:-1] + xedges[1:]) / 2
    y = (yedges[:-1] + yedges[1:]) / 2
    contourf(x, y, H.T, levels)
