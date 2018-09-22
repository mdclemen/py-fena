# Example 6.16
# Two-Dimensional poisson equation
from pylab import*
from mpl_toolkits.mplot3d import Axes3D
import time
# Initialization
clf()

dat = load("ch6ex16_491.npz") # load 37 element mesh from stored data
nodes = dat["nodes"]
elements = dat["elements"]
nin = dat["nin"] # no. of internal nodes
nex = dat["nex"] # no. of external nodes
nce = dat["nce"] # no. of nodes on chamfer (with nonzero bcs)
nn = nodes.shape[0]
mm, ln = elements.shape # 3 nodes per element
t0 = time.time()

As = zeros(mm) # store area of each element/cell
for m in range(mm):
    # each element/cell has coordinates for each node, expressed as i,j,k in the example, so use 0,1,2
    xj = nodes[elements[m,1],0]
    xi = nodes[elements[m,0],0]
    yk = nodes[elements[m,2],1]
    yi = nodes[elements[m,0],1]
    yj = nodes[elements[m,1],1]
    xk = nodes[elements[m,2],0]
    As[m] = abs((xj - xi)*(yk - yi) - (yj - yi)*(xk - xi))/2.

# construct basis functions
def bf(xi, xj, xk, yi, yj, yk):
    den = ((xi - xj)*(yk - yj) - (yi - yj)*(xk - xj))
    return array([(yk - yj)/den, -(xk - xj)/den])

phis = zeros([mm,nn,2]) # order is element/cell, node, [a, b] (basis function coefficients, 6.95)
for m in range(mm): # loop through every element/cell
    for i in elements[m,:]: # loop through every node that composes the element/cell
        # idx is changed for each node that composes a cell/element, in order for that node to be indexed first (i.e. switch i,j,k accordingly). See eqn 6.94
        if i == elements[m,0]:
            idx = [0, 1, 2]
        elif i == elements[m,1]:
            idx = [1, 2, 0]
        elif i == elements[m,2]:
            idx = [2, 0, 1]
        xi = nodes[elements[m,idx]][0,0]
        xj = nodes[elements[m,idx]][1,0]
        xk = nodes[elements[m,idx]][2,0]
        yi = nodes[elements[m,idx]][0,1]
        yj = nodes[elements[m,idx]][1,1]
        yk = nodes[elements[m,idx]][2,1]
        phis[m,i,:] = bf(xi, xj, xk, yi, yj, yk)

K = zeros([nn,nn]) # "stiffness" matrix, 6.92 & 6.96
for i in range(nn):
    for j in range(nn):
        for m in range(mm):
            if i in elements[m,:] and j in elements[m,:]:
                K[i,j] += As[m]*(phis[m,i,0]*phis[m,j,0] + phis[m,i,1]*phis[m,j,1]) # Sum(A*(ami*amj + bmi*bmj))

M = zeros([nn,nn]) # "mass" matrix, 6.98
for i in range(nn):
    for j in range(nn):
        val = 0.
        for m in range(mm):
            if i in elements[m,:] and j in elements[m,:]:
                val += As[m]
        if i == j:
            M[i,j] = (1./6.)*val
        else:
            M[i,j] = (1./12.)*val

# boundary condition on chamfered section and analytical solution
bc = lambda x, y: sin(pi*x)*sin(pi*y)

# q function
qfun = lambda x, y: -2.*(pi**2)*sin(pi*x)*sin(pi*y)
q = qfun(nodes[:,0], nodes[:,1])

ui = hstack((zeros(nex - (nn - nce)), bc(nodes[nce::,0], nodes[nce::,1])))
print("solving for u...")
u = linalg.solve(-K[0:nin,0:nin], dot(K[0:nin,nin::], ui) + dot(M[0:nin,:], q))
u = hstack((u, ui))
tf = time.time()
print("Solution took: %.2f mins, for %d nodes, %d elements" % ((tf - t0)/60., nn, mm))

np.savez("ch6ex16_soln_%de" % mm, usoln = u)
