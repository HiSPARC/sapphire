import tables
from scipy import optimize

from landau import Scintillator

try:
    data
except NameError:
    data = tables.openFile('kascade.h5', 'r')

ph0 = data.root.coincidences.events[:]['pulseheights'][:,0]
s = Scintillator()

n, bins, patches = hist(ph0, bins=linspace(0, 2000, 101), histtype='step')
nx = bins[:-1] + .5 * (bins[1] - bins[0])
x = linspace(-2000, 2000, 200)
y = interp(x, nx, n)
p = optimize.fmin(s.residuals, (3.38 / 380., 10 ** 4, 1), (x, y, 250, 500))
plot(x, s.conv_landau(x, *p))

f = lambda x, N, a: N * x ** -a
x2 = x.compress((0 <= x) & (x < 100))
y2 = y.compress((0 <= x) & (x < 100))
popt, pcov = optimize.curve_fit(f, x2, y2, sigma=y2)
plot(x2, f(x2, *popt))

xscale('log')
yscale('log')
xlim(10, 2 * 10 ** 3)
ylim(ymin=100)
