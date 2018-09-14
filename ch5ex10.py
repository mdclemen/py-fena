# Example 5.10
# V Cycle Multigrid for the Poisson Equation
from pylab import*
from mgv import*
# Setup
clf()

N = 32
n = N + 1

x0 = -1.
xN = 1.
y0 = -1.
yN = 1.

dx = (xN-x0)/N
dy = (yN-y0)/N

x = linspace(x0, xN, N+1)
y = linspace(y0, yN, N+1)

X, Y = meshgrid(x,y)

Q = 2.*(2. - np.power(X, 2) - np.power(Y, 2))

# True solution
Ut = multiply((np.power(X, 2) - 1.), (np.power(Y, 2) - 1.))

# Multi-grid V-cycles
U = asmatrix(zeros([n,n]))
B = -Q

num_vcycles = 4
err = []

for v in range(num_vcycles):
    U = mgv2d(U,B,3)
    err.append(abs(U-Ut).max() / Ut.max())

# Plot
figure(1)
semilogy(err,"bo-", mfc = "none")
xlabel("iteration")
ylabel("error")
xlim([0., num_vcycles])
show()
