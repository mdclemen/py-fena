# plot Example 6.17, plot in a separate routine after solutions are computed
# ch6ex17.py must be run with ch6ex17_1373.npz prior to running this script
from pylab import*
clf()

u = load("ch6ex17_soln_2544e.npz")["usoln"]
dat = load("ch6ex17_1373.npz")
nodes = dat["nodes"]
elements = dat["elements"]
m, nm = elements.shape
no = dat["no"]
ni = dat["ni"]

# plot domain mesh
subplot(121)
ax = gca()

for i in range(m):
    x = zeros(nm+1)
    y = zeros(nm+1)

    for j in range(nm):
        x[j] = nodes[elements[i,j], 0]
        y[j] = nodes[elements[i,j], 1]
    x[nm] = nodes[elements[i,0],0] # restart at the begining of cell, so a triangle is plotted
    y[nm] = nodes[elements[i,0],1]
    for j in range(nm):
        ax.plot(x[j:j+2], y[j:j+2], c = "k", lw = .5, zorder = 0)
ax.set_xticks([-1., -.5, 0., .5, 1.])
ax.set_yticks([-1., -.6, -.2, .2, .6, 1.])
ax.set_aspect("equal")

subplot(122)
ax = gca()
# plot boundaries
ax.plot(nodes[0:no,0], nodes[0:no,1], "k")
ax.plot(nodes[no:no+ni,0], nodes[no:no+ni,1], "k")
# plot solution
lvls = [.1, .2, .3, .4, .5, .6, .7, .8, .9]
mlocs = [(.55, -.7), (.5, -.5), (.4, -.35), (.4, -.15)]
CS = ax.tricontour(nodes[:,0], nodes[:,1], u, colors = "k", linewidths = .7, triangles = elements, levels = lvls)
clabel(CS, inline = 1, fontsize = 8, fmt = "%.1f", manual = mlocs)
ax.set_xticks([-1., -.5, 0., .5, 1.])
ax.set_yticks([-1., -.6, -.2, .2, .6, 1.])
ax.set_aspect("equal")
show()
