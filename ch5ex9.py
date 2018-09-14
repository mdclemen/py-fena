# Example 5.9
# One-Dimensional V Cycle Multigrid
#
from pylab import*
from mgv import*
# Setup
clf()

N = 64
n = N+1
x0 = 0
xN = 1
h = 1./N

X = asmatrix(linspace(x0, xN, N+1)).T
b = (sin(pi*X) + sin(16*pi*X))/2.

# V-cycles
u = asmatrix(zeros([N+1,1]))
num_vcycles = 16
R = zeros(num_vcycles)

for v in range(num_vcycles):
    u = mgv1n(u,b)
    r = asmatrix(zeros([n,1]))
    r[1:n-1,0] = b[1:n-1,0] + (2.*u[1:n-1,0]-u[0:n-2,0]-u[2:n,0]) / h**2
    R[v] = max(abs(r))

R1 = copy(R)
nv1 = num_vcycles

# V-cycles - simple restriction
u = asmatrix(zeros([N+1,1]))
num_vcycles = 28
R = zeros(num_vcycles)

for v in range(num_vcycles):
    u = mgv1r(u,b)
    r = zeros([n,1])
    r[1:n-1] = b[1:n-1] + (2.*u[1:n-1]-u[0:n-2]-u[2:n]) / h**2
    R[v] = max(abs(r))

# Plot
figure(1)
semilogy(range(nv1), R1, "bo-", label = "Average Restriction", mfc = "none")
semilogy(range(num_vcycles), R, "ro-", label = "Simple Restriction", mfc = "none")

xlabel("Iteration")
ylabel("max |r_i|")
legend(loc = 0)
grid("on")
show()
