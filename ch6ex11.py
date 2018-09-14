# Example 6.11
# Calculation of Integrals Using Discrete Chebyshev Transform
from pylab import*
from dtrans import*
# Initialization

clf()
N  = 16;
example = "b"

# Functions we integrate
if example == "a":
    u  = lambda y: (pi - 1.)/4. * sin(0.5*(pi - 1.)*y + 0.5*(pi + 1.))/np.power(0.5*(pi-1.)*y + 0.5*(pi + 1.), 3)

elif example == "b":
    u  = lambda y: 7./2. * log(3.5*y + 4.5)/(3.5*y + 4.5)

# Get the Chebyshev Transform
j  = r_[N:-1:-1]
X  = cos(j*pi/N)
U  = u(X)

a  = dct1(asmatrix(U).T)
a  = vstack((a, 0.))
d  = zeros([N+2,1])

# Numerical Integration
for j in range(2,N+1):
    d[j,0] = (a[j-1,0] - a[j+1,0])/(2.*(j))
d[1,0]   = (2.*a[0,0] - a[2,0])/2.
d[N+1] = a[N,0]/(N+1.)

for j in range(1,N+2):
    d[0,0] = d[0,0] + (-1.)**(j+1) * d[j,0]

# I(x) is know directly in x = 1

print("I = %.15f" % d.sum())
