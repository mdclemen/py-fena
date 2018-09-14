# Example 1.1
# Lagrange Polynomial Interpolation
from pylab import*
from interps import*
# Define vectors

# X : interpolation points
# Y : value of f(X)
# x : points where we want an evaluation of P(x),
#     where P is the interpolator polynomial
x = r_[-1.:1. + .01:.01]
X = r_[-1.:1. + .20:.20]

y = 1./(1. + 25.*np.power(x, 2))
Y = 1./(1. + 25.*np.power(X, 2))

# Lagrange Interpolation
pol = lagrange_interp(X, Y, x)

# Generate plots
plot(x, pol, "k", label = "Lagrange Polynomial")
plot(x, y, "k--", label = "Expected behavior")
plot(X, Y, "k.", label = "Data Points")
legend(loc = 0)
grid("on")

# Evaluate P(.7)
print("p7 =", lagrange_interp(X,Y,.7))
show()
