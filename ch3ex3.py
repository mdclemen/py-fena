# Example 3.3
# Adaptive quadrature
from pylab import*
from scipy.integrate import quad
fun = lambda x: 10.*exp(-50.*abs(x)) - 0.01/(np.power(x - 0.5, 2) + 0.001) + 5.*sin(5.*x)

tol = 1.e-10 #relative tolerance
I   = quad(fun,-1.00, 1.00, epsrel = tol, full_output = 1)
print("I = %.15f" % I[0])
