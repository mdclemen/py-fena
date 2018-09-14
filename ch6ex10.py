# Example 6.10
# Calculation of derivatives using discrete Chebychev Transform
from pylab import*
from dtrans import*
# Initialization
clf()

example = "b"
N = 5

if example == "a":
        u  = lambda x: np.power(x, 4)
        du = lambda x:  4.*np.power(x, 3)
        
elif example == "b":
        u  = lambda x: 4.*multiply(np.power(x, 2), multiply((1. - np.power(x, 2)), exp(-x/2.)))
        du = lambda x: multiply(8.*x - 2.*np.power(x, 2) - 16.*np.power(x, 3) + 2.*np.power(x, 4), exp(-x/2.))

# Using Chebyshev Polynomial :

j  = r_[N:-1:-1]
X1 = cos(j*pi/N)
U1 = u(X1)

a  = dct1(asmatrix(U1).T)
b  = asmatrix(zeros([N+2,1]))

for j in range(N-1,0,-1):
    b[j,0] = b[j+2,0] + 2.*(j+1)*a[j+1,0]

b[0,0] = a[1,0] + 0.5*b[2,0]

b  = -b[0:N+1,0]
b[0,0] = b[0,0]*2.

Ud1 = dct1(b)*N   
Ud1[1:N,0] = Ud1[1:N,0]/2.

# Using Finite differences :

dX = 2./N
X2 = r_[-1.:1. + dX:dX]
U2 = u(X2)

A  = -diag(ones(N), -1) + diag(ones(N), 1)
A[0,0:3] = array([-3., 4., -1.])
A[N,N-2:N+1] = array([1., -4., 3.])
A  = A/(2.*dX)

Ud2 = dot(A, U2)

x = r_[-1.:1. + .01:.01]
plot(X1, Ud1, "ks", mfc = "none", label = "Chebyshev, N = %d" % N)
plot(X2, Ud2, "k*", mfc = "none", label = "Central FD, N = %d" % N)
plot(x, du(x),"k-", label = "Exact")
legend(loc = 0)
xlabel("x"); ylabel("u`")
grid("on")
show()
