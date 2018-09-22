# plot Example 6.16, plot in a separate routine after solutions are computed
# ch6ex16.py must be run with ch6ex16_37.npz and ch6ex16_491.npz prior to running this script
from pylab import*
clf()

u0 = load("ch6ex16_soln_37e.npz")["usoln"]
dat = load("ch6ex16_37.npz")
nodes0 = dat["nodes"]
elements0 = dat["elements"]
m0 = elements0.shape[0]

# domain vertices
domx = array([0., 1., 1., .5, 0., 0.])
domy = array([0., 0., .5, 1., 1., 0.])

# plot 37 element solution
subplot(121)
ax = gca()
ax.plot(domx, domy, "k", lw = .8)
lvls = [.1, .2, .3, .4, .5, .6, .7, .8]
mlocs = [(.15, .05), (.25, .1), (.32, .15), (.38, .2), (.4, .22), (.46, .27), (.50, .3), (.48, .42)]
CS = ax.tricontour(nodes0[:,0], nodes0[:,1], u0, colors = "k", linewidths = .7, triangles = elements0, levels = lvls)
clabel(CS, inline = 1, fontsize = 8, fmt = "%.1f", manual = mlocs)
ax.set_aspect("equal")
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("%d element solution" % m0, fontsize = 9)

u1 = load("ch6ex16_soln_491e.npz")["usoln"]
dat = load("ch6ex16_491.npz")
nodes1 = dat["nodes"]
elements1 = dat["elements"]
m1 = elements1.shape[0]
# generate analytic solution
soln = lambda x, y: sin(pi*x)*sin(pi*y)
x = linspace(0., 1., 100)
y = linspace(0., 1., 100)
X, Y = meshgrid(x, y)
U = soln(X, Y)
fun = lambda x: -x + 1.5
for i in range(x.shape[0]):
    for j in range(y.shape[0]):
        if X[0,i] > .5:
            ymax = fun(X[0,i])
            if Y[j,0] >= ymax:
                U[i,j] = 0.

# plot 37/491 elements soln and analytic soln
subplot(122)
ax = gca()
ax.plot(domx, domy, "k", lw = .8)
mlocs = [(.5, .1), (.5, .24)]
lvls = [.2, .7]
CS0 = ax.tricontour(nodes0[:,0], nodes0[:,1], u0, colors = "k", linewidths = .7, triangles = elements0, levels = lvls, linestyles = "--")
CS0.collections[0].set_label("%d elements" % m0)
CS1 = ax.tricontour(nodes1[:,0], nodes1[:,1], u1, colors = "k", linewidths = .7, triangles = elements1, levels = lvls, linestyles = "dotted")
CS1.collections[0].set_label("%d elements" % m1)
CSa = ax.contour(X, Y, U, levels = lvls, colors = "k", linewidths = .7)
clabel(CSa, inline = 1, fontsize = 8, manual = mlocs, label = "exact")
CSa.collections[0].set_label("exact")
ax.set_aspect("equal")
ax.set_xticks([])
ax.set_yticks([])
ax.legend(loc = "upper right", fontsize = 8)
show()
