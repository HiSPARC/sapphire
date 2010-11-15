import tables
from kascade_coincidences import *

from pylab import *


if __name__ == '__main__':
    try:
        data
        h
        k
    except NameError:
        data = tables.openFile('kascade.h5', 'r')
        h = array(data.root.datasets.h.read())
        k = array(data.root.datasets.knew.read())

    c = array(search_coincidences(h, k, timeshift=-13.180220188))
    cr = array(search_coincidences_right(h, k, timeshift=-13.180220188))
