# Example 6.9
# Calculation of the discrete Chebyshev coefficients
from pylab import*
from dtrans import dct1
# Initialization
clf()

N = 8
j = r_[N:-1:-1]
X = asmatrix(cos(j*pi/N)).T

# Coefficients of x^4

print(dct1(np.power(X, 4)))

# u(x) = 4*(x²-x^4)*e(-x²/2)

u = lambda x: 4.*multiply((np.power(x, 2)), multiply((1. - np.power(x, 2)), exp(-x/2.)))
U = u(X)

plot(r_[0:N+1], abs(dct1(U)), "k-o", mfc = "none", lw = 2.)
xlabel("n"); ylabel("|a_n|")
grid("on")
show()
