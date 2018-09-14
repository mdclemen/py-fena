from pylab import*

def tri_diag(A, rhs):
# x = tri_diag(A,b) takes a tri-diagonal square matrix A and a vector b and 
# returns the solution of Ax=b.

    n, n = shape(A)
    a = zeros(n)
    b = zeros(n)
    c = zeros(n)
    alpha = zeros(n)
    beta = zeros(n)
    y = zeros(n)
    x = zeros(n)

    for i in range(n):
        b[i] = A[i,i]
    for i in range(n-1):
        c[i] = A[i,i+1]
    for i in range(1, n):
        a[i] = A[i,i-1]

    # solve for the entries of L and U so that LU = A:
    beta[0] = b[0]
    for j in range(1, n):
        alpha[j] = a[j]/beta[j-1]
        beta[j] = b[j] - alpha[j]*c[j-1]

    # solve Ly = b
    y[0] = rhs[0]
    for j in range(1, n):
        y[j] = rhs[j] - alpha[j]*y[j-1]

    # solve Ux = y
    x[n-1] = y[n-1]/beta[n-1]
    for j in range(1,n):
        x[n-j-1] = (y[n-j-1]-c[n-j-1]*x[n-j])/beta[n-j-1]

    return x

def makeAb(Nelm, ffunc):
# Nelm = number of elements
# ffunc = string with f function; 

    Delta=1./Nelm
    N = Nelm-1
    A = zeros([N,N])
    b = zeros(N)
    for j in range(N):
        A[j,j] = -(2./Delta) + 2./3.*Delta
    for j in range(N-1):
        A[j,j+1] = 1./Delta + Delta/6.
        A[j+1,j] = 1./Delta + Delta/6.

    x = linspace(0.,1. ,Nelm+1)
    b[0] = Delta/6.*ffunc(x[0]) + 2./3.*Delta*ffunc(x[1]) + Delta/6.*ffunc(x[2])
    b[N-1] = Delta/6.*ffunc(x[N-1]) + 2./3.*Delta*ffunc(x[N]) + Delta/6.*ffunc(x[N+1])
    for j in range(1, N-1):
        b[j] = Delta/6.*ffunc(x[j]) + 2./3.*Delta*ffunc(x[j+1]) + Delta/6.*ffunc(x[j+2])

    xgrid=x[1:Nelm]

    return A, b, xgrid

def makef(phi, beta):
# [A, f]=makef(phi, beta);

    N = size(phi)
    f = zeros(N)
    A = zeros([N,N])

    # Homogeneous boundary conditions
    phi0 = 0.
    phiN1 = 0.

    # Construct the vector f
    f[0] = (1./6. + beta)*phi0 + (2./3. - 2.*beta)*phi[0] + (1./6. + beta)*phi[1]
    f[N-1] = (1./6. + beta)*phi[N-2] + (2./3. - 2.*beta)*phi[N-1] + (1./6. + beta)*phiN1
    for i in range(1, N-1):
        f[i] = (1./6. + beta)*phi[i-1] + (2./3. - 2.*beta)*phi[i] + (1./6. + beta)*phi[i+1]

    # Construct the matrix A
    for j in range(N):
        A[j,j] = (2./3. + 2.*beta)
    for j in range(N-1):
        A[j,j+1] = 1./6. - beta
        A[j+1,j] = 1./6. - beta

    return A, f
