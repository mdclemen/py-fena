#20180921 MDC, reorder nodes from fenics gmsh out so that internal nodes are 0:n and external nodes are n:end, with the chamfered nodes at the very end
from pylab import*

dat = load("ch6ex16_37.npz.orig")
nodes = dat["nodes"]
elements = dat["cells"]
nn = nodes.shape[0]
ne = elements.shape[0]

# remove redundant nodes that are not in any elements
realNodes = []
fakeNodes = []
nodeMap = {}
nodeMap[elements.max()] = 0
for i in range(nn):
    if i in elements:
        realNodes.append(i)
    else:
        fakeNodes.append(i)
    nodeMap[i] = i
print("removing %d redundant nodes" % len(fakeNodes))
# remap nodal indices
rn = len(realNodes)
newNodes = zeros([rn, nodes.shape[1]])
for i in range(rn):
    newNodes[i,:] = nodes[realNodes[i],:]
    nodeMap[realNodes[i]] = i
newElements = zeros_like(elements)
for i in range(ne):
    for j in range(elements.shape[1]):
        newElements[i,j] = nodeMap[elements[i,j]]
nodes = newNodes # overwrite old nodes and elements
elements = newElements
nn = nodes.shape[0]

nex = 0 # number of exterior nodes
nc = 0 # number of nodes on chamfer
newNodes = zeros_like(nodes)
nodeMap = {}
nodeMap[elements.max()] = 0

fun = lambda x: -x + 1.5

for i in range(nn):
    # identify boundary points
    if allclose(nodes[i,1], fun(nodes[i,0]), atol = 1.e-2): # if on chamfered boundary
        nex += 1
        nc += 1
    elif nodes[i,0] == 0. or nodes[i,0] == 1. or nodes[i,1] == 0. or nodes[i,1] == 1.:
        nex += 1
    nodeMap[i] = i
nin = nn - nex #number of interior nodes
nce = nn - nc
ii = 0
jj = 0
ic = 0
for i in range(nn):
    if allclose(nodes[i,1], fun(nodes[i,0]), atol = 1.e-2):
        newNodes[nce+ic,:] = nodes[i,:]
        nodeMap[i] = nce+ic
        ic += 1
    elif nodes[i,0] == 0. or nodes[i,0] == 1. or nodes[i,1] == 0. or nodes[i,1] == 1.:
        newNodes[nin+ii,:] = nodes[i,:]
        nodeMap[i] = nin+ii
        ii += 1
    else:
        newNodes[jj,:] = nodes[i,:]
        nodeMap[i] = jj
        jj += 1
newElements = zeros_like(elements)
for i in range(ne):
    for j in range(elements.shape[1]):
        newElements[i,j] = nodeMap[elements[i,j]]

savez("ch6ex16_37", nodes = newNodes, elements = newElements, nex = nex, nin = nin, nce = nce)
