# Example 5.6
# Approximate Factorization for the Heat Equation
#
# This example uses Scipy's sparse matrices as the operators.  This makes
# the code short and clean.  But, care must be taken in applying the
# correct operator to the correct data.  Matlab's meshgrid function stores
# x-values horizontally and y-values vertically in the data.  When applying
# a matrix operator to data you are operating on the columns.
#
# For example, say T is the data array.  Ay is an operator on y-data and Ax
# is an operator on x-data.  The application must be done in this way:
#   Ay*T  - to operate on the y-data
#   Ax*T' - data must be transposed before applying the x operator
#
# Note 1: An alternative to this method would be to store all of the data
# in a 1-dimensional array.  However, this complicates the construction of
# the operators and plotting in matlab.
#
# Note 2: Comments in this code use the following terms
#  y-major: when y-data is organized in columns, the defualt state
#  x-major: when x-data is organized in columns, usually after a transpose
from pylab import*
from scipy.sparse import diags, eye as speye
from mpl_toolkits.mplot3d import Axes3D
# Set up
clf()
M = 20  # points in the x direction
N = 20  # points in the y direction

x0 = -1.
xM = 1.
y0 = -1.
yN = 1.

dx = (xM-x0)/M
dy = (yN-y0)/N

x = asmatrix(linspace(x0, xM, M+1)).T
y = asmatrix(linspace(y0, yN, N+1)).T

X, Y = meshgrid(x, y)

q = 2.*(2. - np.power(X, 2) - np.power(Y, 2))

# operators
Ax = 1./(dx**2)*asmatrix(diags(array([ones(N-2), -2.*ones(N-1), ones(N-2)]), [-1, 0, 1]).toarray())
Ay = 1./(dy**2)*asmatrix(diags(array([ones(M-2), -2.*ones(M-1), ones(M-2)]), [-1, 0, 1]).toarray())

Ix = speye(M-1)
Iy = speye(N-1)

# Time advancement 1
dt = 0.05             # time step
T = asmatrix(zeros([N+1,M+1]))    # initial condition

t_final = 1.0         # final time
time = r_[0.:t_final+dt:dt]   # time array
pt = array([0.0, 0.25, 1.0]) # desired plot times
pn = len(pt)       # number of desired plots
pc = 0             # plot counter
rt = zeros([1,pn])      # actual plot times
S = zeros([N+1,M+1,pn])  # solution storage

for t in time:
    # plot storage
    if (t >= pt[pc]):
        S[:,:,pc] = T
        rt[0,pc] = t
        pc += 1
        if (pc > pn):
            break
    
    r1 = (Iy + dt/2.*Ay) * T[1:N,1:M];  # r1 is y-major
    r2 = (Ix + dt/2.*Ax) * r1.T         # r2 is x-major
    r = r2 + dt*q[1:N,1:M].T            # r is x-major
    
    T1 = linalg.solve(Ix - dt/2.*Ax, r)           # T1 is x-major
    T[1:N,1:M] = linalg.solve(Iy - dt/2.*Ay, T1.T)# T is y-major

# Plot
fig = figure(1)
ax = fig.add_subplot(111, projection = "3d")
ax.plot_surface(X, Y, S[:,:,0], cmap = "viridis")
ax.set_zlim([0., 1.])
title("t = %.2f" % rt[0,0])

fig = figure(2)
ax = fig.add_subplot(111, projection = "3d")
ax.plot_surface(X, Y, S[:,:,1], cmap = "viridis")
ax.set_zlim([0., 1.])
title("t = %.2f" % rt[0,1])

fig = figure(3)
ax = fig.add_subplot(111, projection = "3d")
ax.plot_surface(X, Y, S[:,:,2], cmap = "viridis")
ax.set_zlim([0., 1.])
title("t = %.2f" % rt[0,2])

# Compute error

# true solution
P = multiply((np.power(X, 2) - 1.), (np.power(Y, 2) - 1.))

# max pointwise difference
mpd1 = (abs(P-S[:,:,2])).max()
print("mpd1 = %e" % mpd1)
# max pointwise percentage
mpp1 = (abs(P[1:N,1:M]-S[1:N,1:M,2])/P[1:N,1:M]).max()
print("mpp1 = %e" % mpp1)
# percentage error volume
pev1 = (abs(P-S[:,:,2])).sum()/P.sum()
print("pev1 = %e" % pev1)

# Time advancement 2
dt = 1.               # time step
T = asmatrix(zeros([N+1,M+1]))    # initial condition

t_final = 5.0         # final time
time = r_[0.:t_final+dt:dt]   # time array
pt = [5.0]            # desired plot times
pn = len(pt)          # number of desired plots
pc = 0                # plot counter
rt = zeros([1,pn])      # actual plot times
S = zeros([N+1,M+1,pn]) # solution storage

for t in time:
    # plot storage
    if (t >= pt[pc]):
        S[:,:,pc] = T
        rt[0,pc] = t
        pc += 1
        if (pc > pn):
            break
    
    r1 = (Iy + dt/2.*Ay) * T[1:N,1:M]  # r1 is y-major
    r2 = (Ix + dt/2.*Ax) * r1.T        # r2 is x-major
    r = r2 + dt*q[1:N,1:M].T           # r is x-major
    
    T1 = linalg.solve(Ix - dt/2*Ax, r) # T1 is x-major
    T[1:N,1:M] = linalg.solve(Iy - dt/2*Ay, T1.T) # T is y-major

# Compute error

# true solution
P = multiply((np.power(X, 2) - 1.), (np.power(Y, 2) - 1.))

# max pointwise difference
mpd2 = (abs(P-S[:,:,0])).max()
print("mpd2 = %e" % mpd2)
# max pointwise percentage
mpp2 = (abs(P[1:N,1:M]-S[1:N,1:M,0])/P[1:N,1:M]).max()
print("mpp2 = %e" % mpp2)
# percentage error volume
pev2 = (abs(P-S[:,:,0])).sum()/P.sum()
print("pev2 = %e" % pev2)
show()
