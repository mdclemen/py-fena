# Example 3.4
# Gauss - Legendre quadrature
from pylab import*
from gauss import*
a   = 1.
b   = 8.
fun = lambda x: log(x)/x

n   = 5
x, w = gauss_leg(a,b,n);
print("I = %.15f" % dot(w, fun(x)))
