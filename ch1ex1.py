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

# Using only numpy
def prod(x, xs, nn, idx):
    res = 1.
    for i in range(nn):
        if i != idx:
            res *= (x - xs[i])/(xs[idx] - xs[i])
    return res

n, = X.shape
m, = x.shape
yy = zeros(m)
for i in range(m):
    for j in range(n):
        yy[i] += Y[j]*prod(x[i], X, n, j) # Moin eq. 1.2


# Generate plots
plot(x, pol, "k", label = "Lagrange Polynomial")
plot(x, yy, "k-.", label = "Lagrange Polynomial (2nd method)")
plot(x, y, "k--", label = "Expected behavior")
plot(X, Y, "k.", label = "Data Points")
legend(loc = 0)
grid("on")

# Evaluate P(.7)
print("p7 =", lagrange_interp(X,Y,.7))
show()
