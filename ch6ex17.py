# Example 6.17
# Steady heat equation
from pylab import*
from mpl_toolkits.mplot3d import Axes3D
import time, argparse, subprocess
parser = argparse.ArgumentParser(description = "Solve poisson's equation on unit square.")
parser.add_argument("--useC", dest = "useC", action = "store_true", default = False, help = "Use C functions to compute K and M.")
parser.add_argument("--useOCL", dest = "useOCL", action = "store_true", default = False, help = "Use OpenCL functions to compute K and M.")
args, unknown = parser.parse_known_args()
# Initialization
clf()

dat = load("ch6ex17_1373.npz") # load element mesh from stored data
nodes = dat["nodes"]
elements = dat["elements"]
nn = nodes.shape[0]
mm, ln = elements.shape # 3 nodes per element
no = dat["no"] # number of outer boundary nodes
ni = dat["ni"] # number of inner boundary nodes
nex = no + ni # no of external/bc nodes
nin = nn - nex # no of internal/unknown nodes
if args.useC:
    from ctypes import c_int, c_void_p, CDLL
    try:
        funcs = CDLL("./funcs.so")
    except:
        subprocess.call(["make", "all"])
        funcs = CDLL("./funcs.so")        
    makeK = funcs.makeK
    makeM = funcs.makeM
elif args.useOCL:
    import pyopencl as cl
    ctx = cl.create_some_context()
    queue = cl.CommandQueue(ctx)
    mf = cl.mem_flags
    krnl_src = open("funcs.cl", "r").read()
    funcs = cl.Program(ctx, krnl_src).build("-I./")
print("Calculating areas...")
t0 = time.time() # start timing
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
t1 = time.time()
print("Area computation: %.2f mins" % ((t1 - t0)/60.))
# construct basis functions
def bf(xi, xj, xk, yi, yj, yk):
    den = ((xi - xj)*(yk - yj) - (yi - yj)*(xk - xj))
    return array([(yk - yj)/den, -(xk - xj)/den])
print("Calculating basis functions...")
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
t2 = time.time()
print("Basis function computation: %.2f mins" % ((t2 - t1)/60.))
print("Calculating stiffnes matrix K...")
K = zeros([nn,nn]) # "stiffness" matrix, 6.92 & 6.96
if args.useC:
    makeK(nn, mm, c_void_p(K.ctypes.data), c_void_p(As.ctypes.data), c_void_p(phis.ctypes.data))
elif args.useOCL:
    d_K = cl.Buffer(ctx, mf.WRITE_ONLY | mf.COPY_HOST_PTR, hostbuf = K)
    d_As = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf = As)
    d_phis = cl.Buffer(ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf = phis)
    funcs.makeK(queue, (nn,), None, int32(nn), int32(mm), d_K, d_As, d_phis).wait()
    cl.enqueue_copy(queue, K, d_K)
else:
    for i in range(nn):
        for j in range(nn):
            for m in range(mm):
                # if i in elements[m,:] and j in elements[m,:]:
                K[i,j] += As[m]*((phis[m,i,0]*phis[m,j,0]) + (phis[m,i,1]*phis[m,j,1])) # Sum(A*(ami*amj + bmi*bmj))
t3 = time.time()
print("Stiffness matrix computation: %.2f mins" % ((t3 - t2)/60.))
# q function is 0, so calculation of Mass matrix is not necessary

ui = hstack((zeros(no), ones(ni)))# outer boundary node bcs, inner boundary node bcs
print("Solving for u...")
u = linalg.solve(-K[nex:nn,nex:nn], dot(K[nex:nn,0:nex], ui))
u = hstack((ui, u))
tf = time.time()
print("Solution took: %.2f mins, for %d nodes, %d elements" % ((tf - t0)/60., nn, mm))

np.savez("ch6ex17_soln_%de" % mm, usoln = u)
