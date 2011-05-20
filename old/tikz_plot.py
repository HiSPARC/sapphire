from pylab import *

def tikz_2dhist(filename, x, y, bins, use_log=False):
    """Make a 2d histogram and output a tikz plot

    Make a 2d histogram and output a scatter plot where each marker's size
    shows the bin count.  The plot is normed and the size of the markers
    can be chosen to follow log scale.  The output can be included inside
    a tikzpicture environment.  Mind that you need to supply the marker
    style, x and y scales and that you need to manually draw the axis.
    The benefit of not coding that in this function is the ability to
    fine-tune the plot, while still being able to redo an analysis and
    update the plot without having to redo your work.  An example:

        \\begin{tikzpicture}[
            x=8cm/100,y=6cm/10,
            marker/.style={color=gray,line width=0pt,fill},
        ]
        \\input{myplot.tikz}
        \\draw (0,0) rectangle (100,10);
        \\foreach \\x in {0,10,...,100} {
            \\draw (\\x,-2pt) -- (\\x,2pt);
            \\node[below] at (\\x,0) {\\x};
        }
        \\foreach \\y in {0,1,...,10} {
            \\draw (-2pt,\\y) -- (2pt,\\y);
            \\node[left] at (0,\\y) {\\y};
        }
        \\end{tikzpicture}

    """
    H, xedges, yedges = histogram2d(x, y, bins=bins)

    xs = [mean([x1, x2]) for x1, x2 in zip(xedges[:-1], xedges[1:])]
    ys = [mean([y1, y2]) for y1, y2 in zip(yedges[:-1], yedges[1:])]

    xp, yp, sp = [], [], []

    for i in xrange(H.shape[0]):
        for j in xrange(H.shape[1]):
            v = H[i][j]
            if v != 0:
                xp.append(xs[i])
                yp.append(ys[j])
                sp.append(v)

    xr = (xs[1] - xs[0])
    yr = (ys[1] - ys[0])

    x0, x1 = xedges[0], xedges[-1]
    y0, y1 = yedges[0], yedges[-1]

    m = max(sp)
    if use_log:
        sp = [log(x) / log(m) for x in sp]
    else:
        sp = [sqrt(x / m) for x in sp]

    with open(filename, 'w') as file:
        file.write("%% x0, x1: %f, %f\n" % (x0, x1))
        file.write("%% y0, y1: %f, %f\n" % (y0, y1))
        file.write("\\def\\xr{%f}\n" % xr)
        file.write("\\def\\yr{%f}\n" % yr)
        file.write("\\def\\scatterpoint#1#2#3{\n"
                   "  \\draw[marker] (#1-.5*\\xr*#3,#2-.5*\\yr*#3)\n"
                   "    rectangle +(\\xr*#3,\\yr*#3);\n"
                   "}\n")
        for x, y, s in zip(xp, yp, sp):
            file.write("\\scatterpoint{%f}{%f}{%f}\n" % (x, y, s))


if __name__ == '__main__':
    x = linspace(0, 100, 100000)
    y = sqrt(x)

    x += normal(size=x.shape)
    y += normal(size=y.shape)

    tikz_2dhist('myplot.tikz', x, y, bins=(linspace(0, 100, 120),
                linspace(0, 10, 90)), use_log=True)
