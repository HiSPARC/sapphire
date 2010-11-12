import tables
data = tables.openFile('kascade.h5', 'r')
data
print data
h = data.root.datasets.h.read()
k = data.root.datasets.knew.read()
h[:10]
type(h)
h = array(h)
k = array(k)
h[:10,1]
hist(h[:,1], bins=1000)
h[100000:,1] - h[0,1]
h[100000,1] - h[0,1]
(h[100000,1] - h[0,1]) / 1e9
100001 / ((h[100000,1] - h[0,1]) / 1e9)
figure()
hist(k[:,1], bins=1000)
10001 / ((k[10000,1] - k[0,1]) / 1e9)
hr = 100001 / ((h[100000,1] - h[0,1]) / 1e9)
kr = 10001 / ((k[10000,1] - k[0,1]) / 1e9)
hr
kr
hr * kr * 10e-6
import hisparc
help(hisparc.analysis.kascade_coincidences)
c = hisparc.analysis.kascade_coincidences.search_coincidences(h, k, 0, None)
len(c)
c[:10]
type(c)
figure()
hist([x[0] for x in c], bins=200, histtype='step')
hist([x[0]*1e-9 for x in c], bins=200, histtype='step')
clf()
hist([x[0]*1e-9 for x in c], bins=200, histtype='step')
hist([x[0]*1e-9 for x in c], bins=linspace(-1, 1, 100), histtype='step')
xlim(-1, 1)
c[:10]
c2 = hisparc.analysis.kascade_coincidences.search_coincidences(h, k, 0, dtlimit=10000)
len(c2)
c2
from scipy.optimize import curve_fit
hist([x[0]*1e-9 for x in c], bins=linspace(0, 1, 100), histtype='step')
h, bins, p = hist([x[0]*1e-9 for x in c], bins=linspace(0, 1, 100), histtype='step')
xlim(0,1)
ylim(0, 50000)
ylim(0, 20000)
f = lambda x, N, l: N * l * exp(-l * x)
curve_fit(f, bins, h)
len(h)
len(bins)
mx = [(u + v) / 2 for u, v in zip(bins[:-1], bins[1:])]
mx
len(mx)
my = h
len(my)
curve_fit(f, mx, my)
help(curve_fit)
type(mx)
type(my)
mx = array(mx)
curve_fit(f, mx, my)
popt, pcov = curve_fit(f, mx, my)
plot(x, f(x, popt[0], popt[1]))
plot(mx, f(mx, popt[0], popt[1]))
popt
kr
hr
figure()
mx
plot(mx, exp(-mx))
plot(mx, exp(-2 * mx))
plot(mx, exp(-4 * mx))
plot(mx, exp(-6 * mx))
plot(mx, exp(-8 * mx))
