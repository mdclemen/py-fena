# Example 3.5
# Gauss - Hermite quadrature
from pylab import*
from gauss import*
n   = 7
x, w = gauss_her(n)
print("I = %.15f" % dot(w, cos(x)))
