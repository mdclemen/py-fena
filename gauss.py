from pylab import*
from scipy.misc import factorial

def gauss_her(n):
#=====================================================
# Gauss - Hermite quadrature
#=====================================================
# n : number of points for evaluating f
#=====================================================

# Build Hn (Hermite Poly of order n) by recurrence :
    h = zeros([n+1,n+1])
    h[0,0] = 1.
    h[1,0:2] = array([2., 0])
    for k in range(1, n):
        h[k+1,0:k+2] = 2.*hstack((h[k,0:k+1], 0.)) - 2.*(k)*hstack((0., 0., h[k-1,0:k]))

    Hn       = h[n,:]
    Hn_deriv = polyder(Hn)

    x        = roots(Hn)
    w        = 2**(n+1)*factorial(n)*sqrt(pi)/np.power(polyval(Hn_deriv, x), 2)

    return x, w

def gauss_leg(a, b, n):
#=====================================================
# Gauss - Legendre quadrature
#=====================================================
# a : low boundary of integration
# b : high boundary of integration (b>a)
# n : number of points for evaluating f
#=====================================================

    p = zeros([n+1, n+1])
    if (a < b):
        # Build Pn (Legendre Poly of order n) by recurrence :
        p[0,0]   = 1.
        p[1,0:2] = array([1., 0.])
        for k in range(1, n):
            p[k+1,0:k+2] = ((2.*(k + 1) - 1.)*hstack((p[k,0:k+1], 0.)) - (k)*hstack((0., 0., p[k-1,0:k])))/float(k+1.)

        Pn       = p[n,:]
        Pn_deriv = polyder(Pn)

        x        = roots(Pn)
        w        = 2./(multiply(1. - np.power(x, 2), np.power(polyval(Pn_deriv, x), 2)))

        # Go back to interval [a,b] :
        x        = (b + a)/2. + multiply((b - a), x)/2.
        w        = (b - a)/2. * w

        return x, w

    else:
        print("a >= b !!")
        return 0, 0
