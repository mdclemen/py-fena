from pylab import*
from scipy.interpolate import PPoly

def lagrange_interp(X, Y, x):
#=====================================================
# Lagrange Polynomial Interpolation
#=====================================================
# X : interpolation points
# Y : value of f(X)
# x : points where we want an evaluation of P(x),
#     where P is the interpolator polynomial
#=====================================================
    n          = size(X)
    phi        = ones([n, size(x)])
    polynomial = zeros(size(x))

    for i in range(n):
        for j in range(n):
            if (i != j):
               phi[i,:] = multiply(phi[i,:], (x - X[j])/(X[i] - X[j]))

    for i in range(n):
        polynomial = polynomial + Y[i]*phi[i,:]

    return polynomial

def naturalCubicSpline(X, Y, xx):
#===============================================================
# Cubic Spline Interpolation with the Free Run-out End Condition
#===============================================================
# Syntax  
# --------------------------------------------------------------
#   yy = naturalCubicSpline(X,Y,xx)
#   pp = naturalCubicSpline(X,Y)
#
# Input Parameters  
# --------------------------------------------------------------
#   X  : interpolation points
#   Y  : value of f(X)
#   xx : points where we want an evaluation of P(x),
#        where P is the interpolator polynomial
#
# Description
# --------------------------------------------------------------
#   yy = naturalCubicSpline(X,Y,xx) uses a cubic spline interpolation
#   to find yy at the values of the interpolant xx. The end
#   condition is Free run-out (see p.6 in the text). The values 
#   in x must be distinct.
#===============================================================

    # Sort [X,Y] to avoid errors in ppval()
    Xsorted= sort(X, axis = 0)
    I = unravel_index(argsort(X, axis = 0), X.shape)
    Ysorted = Y[I]

    n     = size(Xsorted)
    delta = Xsorted[1:n] - Xsorted[0:n-1]

    # Construct the tri-diagonal matrix to find g'' (with free run-out b.c.)
    matSize = n-2
    M = zeros([matSize, matSize])
    b = zeros(matSize)
    for i in range(matSize):
        if i > 0:
            M[i,i-1] = delta[i]/6.
        M[i,i] = (delta[i] + delta[i+1])/3.
        if i < matSize-1:
            M[i,i+1] = delta[i+1]/6.

        b[i] = (Ysorted[i+2] - Ysorted[i+1])/delta[i+1] - (Ysorted[i+1] - Ysorted[i])/delta[i]

    # Solve the system for g'' and add boundary points
    gpp = linalg.solve(M, b)
    gpp = hstack((0., gpp, 0.))

    # Construct the piecewise polynomials
    coefs = zeros([n-1, 4])
    for i in range(n-1):
        Delta = delta[i]

        coefs[i,0] = (gpp[i+1] - gpp[i])/(6.*Delta)
        coefs[i,1] = gpp[i]/2.
        coefs[i,2] = -(gpp[i+1] + 2.*gpp[i])/6.*Delta + (Ysorted[i+1] - Ysorted[i])/Delta
        coefs[i,3] = Ysorted[i]
    coefs = coefs.T
    pp = PPoly(coefs, Xsorted)

    # Return
    return pp(xx)
