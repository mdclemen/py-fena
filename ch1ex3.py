# Example 1.3
# Interpolation of Runge function using splines
from pylab import*
from interps import*
# Define Vectors
x = r_[-1.:1. + .01:.01]; x[-1] = 1. # stops rounding errors
X = r_[-1.:1. + .20:.2]; X[-1] = 1.

y = 1./(1. + 25.*np.power(x, 2))
Y = 1./(1. + 25.*np.power(X, 2))

# Cubic spline interpolation
# Free end condition
interp = naturalCubicSpline(X, Y, x)

# Using only numpy
n, = X.shape
m, = x.shape
yy = zeros(m)
# Cubic spline interpolation
delta = diff(X) # dim = n-1, uniform here, but doesn't have to be
M = zeros([n-2, n-2])
b = zeros(n-2)
for i in range(n-2): # build tri-diagonal matrix
    if i > 0:
        M[i,i-1] = delta[i]/6.
    M[i,i] = (delta[i] + delta[i+1])/3.
    if i < n-2-1:
        M[i,i+1] = delta[i+1]/6.
    b[i] = (Y[i+2] - Y[i+1])/delta[i+1] - (Y[i+1] - Y[i])/delta[i]

gpp = linalg.solve(M, b)
gpp = hstack([0., gpp, 0.]) # free run-out (natural spline) "end conditions"
for i in range(m): # using Moin equation 1.6
    for j in range(n-1):
        if (X[j] <= x[i] <= X[j+1]): # piecewise
            yy[i] = gpp[j]/6. * (((X[j+1] - x[i])**3)/delta[j] - delta[j]*(X[j+1] - x[i])) +\
                gpp[j+1]/6. * (((x[i] - X[j])**3)/delta[j] - delta[j]*(x[i] - X[j])) +\
                Y[j] * (X[j+1] - x[i])/delta[j] + Y[j+1] * (x[i] - X[j])/delta[j]

# Generates plots
plot(x, interp, "k-", label = "Cubic Spline")
plot(x, yy, "k--", label = "Cubic Spiline (2nd method)")
plot(X, Y, "k.", label = "Data Points")
legend(loc = 0)
xlim([-1., 1.]); ylim([0., 1.2])
grid("on")
show()
