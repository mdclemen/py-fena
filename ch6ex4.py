# Example 6.4
# Poisson Equation with non Homogeneous Boundary Conditions
from pylab import*
from dtrans import dst1
from mpl_toolkits.mplot3d import Axes3D
# Initialization
clf()

M = 32
N = 32
dx= 1./M
dy= 1./N
x = r_[0.:1.+dx:dx]
y = r_[0.:1.+dy:dy]
X,Y = meshgrid(x, y)

Pji = asmatrix(zeros([N+1, M+1]))
Pij = asmatrix(zeros([M+1, N+1]))
Pjk = asmatrix(zeros([M+1, N+1]))

# Define the RHS of modified eq
fun = lambda x, y: asmatrix(30.*(np.power(x, 2) - x) + 30.*(np.power(y, 2) - y) - multiply(4.*pi**2.*(x - 1.), sin(2.*pi*y)))
Qij = fun(X, Y)
Qkj = asmatrix(dst1(Qij.T))
Qjk = asmatrix(Qkj.T)

# Solve using DST
for k in range(M-1):
    #Solve the tridiagonal system for each k
    A = asmatrix(diag(ones(M-2), -1) + diag(ones(M-2), 1) + diag((((dy/dx)**2*(2.*cos(pi*float(k+1.)/float(M)) - 2.)) - 2.)*ones(M-1)))
    b = dy**2 * Qjk[1:M, k+1]
    Pjk[1:M, k+1] = linalg.solve(A, b)

#Inverse transform
Pji = dst1(Pjk.T)*M/2.
Pij = Pji.T

# Solution of the original equation
Sol = Pij - multiply((X - 1.), sin(2.*pi*Y))

# Generates the plots
fig = figure(1)
ax = fig.add_subplot(111, projection = "3d")
ax.plot_surface(X, Y, Sol, cmap = "viridis")
ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("phi")
title("Numerical solution of Poisson equation, ex 6.4")
show()
