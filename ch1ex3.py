# Example 1.3
# Interpolation of Runge function using splines
from pylab import*
from interps import*
# Define Vectors
x = r_[-1.:1. + .01:.01]
X = r_[-1.:1. + .20:.2]

y = 1./(1. + 25.*np.power(x, 2))
Y = 1./(1. + 25.*np.power(X, 2))

# Cubic spline interpolation
# Free end condition
interp = naturalCubicSpline(X, Y, x)

# Generates plots
plot(x, interp, "k-", label = "Cubic Spline")
plot(X, Y, "k.", label = "Data Points")
legend(loc = 0)
xlim([-1., 1.]); ylim([0., 1.2])
grid("on")
show()
