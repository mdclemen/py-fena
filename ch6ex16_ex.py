# Example 6.16
# Two-Dimensional Poisson Equation, using values given in example
from pylab import*

nodes = array([[.3836, .3766], # x1, y1
               [.3736, .6364], # x2, y2
               [.6264, .3966]])# x3, y3

cells = array([[0, 1, 2], # first row of nodal indices
               [1, 2, 0], # [0,2,4] is given but none of those coordinates are given
               [2, 0, 1]])

nn = nodes.shape[0]
mm, ln = cells.shape # 3 nodes per element

# compute area
As = zeros(mm) # store area of each element/cell
for m in range(mm):
    # each element/cell has coordinates for each node, expressed as i,j,k in the example, so use 0,1,2
    xj = nodes[cells[m,1],0]
    xi = nodes[cells[m,0],0]
    yk = nodes[cells[m,2],1]
    yi = nodes[cells[m,0],1]
    yj = nodes[cells[m,1],1]
    xk = nodes[cells[m,2],0]
    As[m] = abs((xj - xi)*(yk - yi) - (yj - yi)*(xk - xi))/2.

print("A_1 = %.5f, 0.03164 given" % As[0])

# construct basis functions
def bf(xi, xj, xk, yi, yj, yk):
    den = ((xi - xj)*(yk - yj) - (yi - yj)*(xk - xj))
    return array([(yk - yj)/den, -(xk - xj)/den])

phis = zeros([mm,nn,2]) # order is element/cell, node, [a, b] (basis function coefficients, 6.95)
for m in range(mm): # loop through every element/cell
    for i in cells[m,:]: # loop through every node that composes the element/cell
        # idx is changed for each node that composes a cell/element, in order for that node to be indexed first (i.e. switch i,j,k accordingly). See eqn 6.94
        if i == cells[m,0]:
            idx = [0, 1, 2]
        elif i == cells[m,1]:
            idx = [1, 2, 0]
        elif i == cells[m,2]:
            idx = [2, 0, 1]
        xi = nodes[cells[m,idx]][0,0]
        xj = nodes[cells[m,idx]][1,0]
        xk = nodes[cells[m,idx]][2,0]
        yi = nodes[cells[m,idx]][0,1]
        yj = nodes[cells[m,idx]][1,1]
        yk = nodes[cells[m,idx]][2,1]
        phis[m,i,:] = bf(xi, xj, xk, yi, yj, yk)

print("(a11, b11) = (%.4f, %.4f), (-3.7987, -3.9950) given\n\
(a12, b12) = (%.4f, %.4f), (-0.3161, 3.8369) given\n" \
% (phis[0,0,0], phis[0,0,1], phis[0,1,0], phis[0,1,1]))
