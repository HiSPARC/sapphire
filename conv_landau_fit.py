from __future__ import division

import tables
from scipy import optimize

from landau import Scintillator

try:
    data
except NameError:
    data = tables.openFile('kascade.h5', 'r')

ph0 = data.root.coincidences.events[:]['pulseheights'][:,0]
s = Scintillator()

figure()
# Fit of convoluted Landau
n, bins, patches = hist(ph0, bins=linspace(0, 2000, 101), histtype='step')
nx = bins[:-1] + .5 * (bins[1] - bins[0])
x = linspace(-2000, 2000, 200)
y = interp(x, nx, n)
p = optimize.fmin(s.residuals, (3.38 / 380., 10 ** 4, 1), (x, y, 250, 500))
plot(x, s.conv_landau(x, *p))

# Fit of gamma spectrum
f = lambda x, N, a: N * x ** -a
x2 = x.compress((0 <= x) & (x < 100))
y2 = y.compress((0 <= x) & (x < 100))
popt, pcov = optimize.curve_fit(f, x2, y2, sigma=y2)
plot(x2, f(x2, *popt))

xscale('log')
yscale('log')
xlim(10, 2 * 10 ** 3)
ylim(ymin=100)
xlabel("Pulseheight [ADC counts]")
ylabel("Counts")
title("Fit of convoluted Landau to data")
savefig("plots/conv_landau_fit.pdf")

figure()
# Charged particle spectrum
step(x, y, where='mid')
yl = s.conv_landau(x, *p)
plot(x, yl)
i = (y <= yl).argmax()
yp = array(yl[:i].tolist() + y[i:].tolist())
step(x, yp, where='mid')
N_T = sum(y.compress(x >= 0))
N_CP = sum(yp.compress(x >= 0))
print "Charged particles: %.2f %% of events" % (N_CP / N_T)

xlim(xmin=0)
yscale('log')
ylim(ymin=10)
ylim(ymin=100)
xlabel("Pulseheight [ADC counts]")
ylabel("Counts")
title("Charged particle part of spectrum")
savefig("plots/conv_landau_charged_particles.pdf")
