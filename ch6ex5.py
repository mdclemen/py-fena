# Example 6.5
# Differentiation Using the Fourier Spectral Method and Second-Order Central Difference Formula
from pylab import*
clf()

# Choice of example a or b : exact solution
example = "b"

if example == "a":
    f  = lambda x: cos(3.*x)
    fd = lambda x: -3.*sin(3.*x)

elif example == "b":
    f  = lambda x: 2.*pi*x - np.power(x, 2)
    fd = lambda x: 2.*pi - 2.*x

# Initialization :
x   = r_[0.:2.*pi + 0.01:0.01]

N   = 8
dX  = 2.*pi/N
X   = r_[0.:2.*pi:dX]
Y   = f(X)

# Fourier Spectral Method :
k   = hstack((r_[0.:N/2.], 0., r_[-N/2. + 1.:0.]))

YY  = fft(Y)/N
YYd = 1j*multiply(k, YY)
Yd1 = ifft(YYd)*N

# 2nd order Finite Difference :
#extension of domain (assuming periodicity)
Y2  = hstack((f(2.*pi-dX), Y, f(2.*pi)))

#differentiation operator for interior points
A   = diag(-ones(N+1), -1) + diag(ones(N+1), 1)
A   = A/(2.*dX)
Yd2 = dot(A, Y2)
Yd2 = Yd2[1:N+1]

# Plot
if example == "a":
    #exact
    plot(x/pi, fd(x), "k-", lw =2., label = "Exact")

    #spectral
    plot(X/pi, Yd1, "ko", label = "Spectral N=%d" % N)

    #FD (for N=16 and N=8)
    plot(X/pi, Yd2,"k:^", lw = 2., label = "FD, N=%d" % N)
    #plot(X/pi, Yd2, "k-.v", lw = 2., label = "FD, N=16")

    xlim([-0.1, 2.2]); ylim([-4., 4.])
        
elif example == "b":
    #exact
    plot(x/pi, fd(x), "k-", lw = 2., label = "Exact")
        
    #spectral
    plot(X/pi, Yd1, "k:o", lw = 2., label = "Spectral N=%d" % N)
    
    #FD
    plot(X/pi, Yd2, "k-.^", lw = 2., label = "FD, N=%d" % N)
    
    xlim([-0.1, 2.0]); ylim([-10., 10.])

legend(loc = 0)
grid("on")
show()
