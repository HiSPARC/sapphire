import tables
import hisparc
from scipy.optimize import curve_fit

try:
    h, k, c, hr, kr
except NameError:
    data = tables.openFile('kascade.h5', 'r')
    h = array(data.root.datasets.h.read())
    k = array(data.root.datasets.knew.read())
    hr = 100001 / ((h[100000,1] - h[0,1]) / 1e9)
    kr = 10001 / ((k[10000,1] - k[0,1]) / 1e9)
    c = hisparc.analysis.kascade_coincidences.search_coincidences(h, k, 0,
                                                                  None)

figure()
counts, bins, p = hist([x[0]*1e-9 for x in c], bins=linspace(0, 1, 100),
                       histtype='step')

f = lambda x, N, l: N * l * exp(-l * x)
mx = array([(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])])
my = counts
popt, pcov = curve_fit(f, mx, my)

plot(mx, f(mx, popt[0], popt[1]))
