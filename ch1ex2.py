# Example 1.2
# Lagrange Polynomial Interpolation of Runge Function
#
# Interpolation of the Runge function using 
# 11 NON equally spaced points on [-1,1]
# (data points are more finely spaced near ends)
from pylab import*
from interps import*
# Define Vectors
x = r_[-1.:1 + .01:.01]
X = array([-1.00, -0.95, -0.81, -0.59, -0.31, 0.00, 0.31, 0.59, 0.81, 0.95, 1.00])

y = 1./(1. + 25.*np.power(x, 2))
Y = 1./(1. + 25.*np.power(X, 2))

# Lagrange Interpolation
pol = lagrange_interp(X, Y, x)

# Generate plots
plot(x, pol, "k", label = "Lagrange Polynomial")
plot(x, y, "k--", label = "Expected behavior")
plot(X, Y, "k.", label = "Data Points")
legend(loc = 0)
xlim([-1., 1.]); ylim([-0.05, 1.2])
grid("on")
show()
