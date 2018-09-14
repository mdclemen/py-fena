# Example 5.8
# Iterative Solution of an Elliptic Equation
from pylab import*
from scipy.sparse import diags, eye as speye
# Setup

N = 20  # points in the each direction

x0 = -1
xN = 1
y0 = -1
yN = 1

dx = (xN-x0)/N
dy = (yN-y0)/N

x = linspace(x0, xN, N+1)
y = linspace(y0, yN, N+1)

X, Y = meshgrid(x,y)

Q = 2.*(2. - np.power(X, 2) - np.power(Y, 2))
Qint = Q[1:N,1:N]  # obtain interior
q = Qint.reshape((N-1)**2,1)        # transform into a vector

# Direct Solve
# construct system.  
# See Applied Numerical Linear Algebra by James Demmel
# footnote on page 275
TN = diags(array([-ones(N-2), 2.*ones(N-1), -ones(N-2)]), [-1,0,1]).toarray()
TNxN = kron(eye(N-1), TN) + kron(TN, eye(N-1))

pd = linalg.solve(1./(dx**2)*TNxN, q)

Pd = zeros([N+1, N+1])
Pd[1:N,1:N] = pd.reshape(N-1,N-1)

# True solution
Pd = multiply(np.power(X, 2) - 1., np.power(Y, 2) - 1.)

tol = 1.e-4

# Point Jacobi
Pj = zeros([N+1,N+1])   # initial guess
Pnew = zeros([N+1,N+1]) # extra storage
err = (abs(Pd[1:N,1:N]-Pj[1:N,1:N])/Pd[1:N,1:N]).max()
niter = 0
Ej = []

while (err > tol):
    niter += 1
    for j in range(1,N):
        for i in range(1,N):
            Pnew[i,j] = (Pj[i+1,j]+Pj[i-1,j]+Pj[i,j+1]+Pj[i,j-1])/4. + dx**2/4.*Q[i,j]
    Pj = copy(Pnew)
    err = (abs(Pd-Pj)).max()/(abs(Pd)).max()
    Ej.append(err)

print("pj_iter = %d" % niter)

# Gauss-Seidel
Pgs = zeros([N+1,N+1])   # initial guess
err = (abs(Pd[1:N,1:N]-Pgs[1:N,1:N])/Pd[1:N,1:N]).max()
niter = 0
Egs = []
while (err > tol):
    niter += 1
    for j in range(1,N):
        for i in range(1,N):
            Pgs[i,j] = (Pgs[i+1,j]+Pgs[i-1,j]+Pgs[i,j+1]+Pgs[i,j-1])/4. + dx**2/4.*Q[i,j]
    err = (abs(Pd-Pgs)).max()/(abs(Pd)).max()
    Egs.append(err)

print("pgs_iter = %d" % niter)

# SOR
Psor = zeros([N+1,N+1])   # initial guess
Pold = zeros([N+1,N+1])   # old solution
omega = 1.8 # SOR parameter
err = (abs(Pd[1:N,1:N]-Psor[1:N,1:N])/Pd[1:N,1:N]).max()
niter = 0
Esor = []
while (err > tol):
    niter += 1
    for j in range(1,N):
        for i in range(1,N):
            Psor[i,j] = omega*((Psor[i+1,j]+Psor[i-1,j]+Psor[i,j+1]+Psor[i,j-1])/4. + dx**2/4.*Q[i,j]) + (1. - omega)*Psor[i,j]
    err = (abs(Pd-Psor)).max()/(abs(Pd)).max()
    Esor.append(err)

print("psor_iter = %d" % niter)
