# Example 6.15
# Unsteady heat equation
from pylab import*
from tri_diag import*

clf()
N=8
alpha = 0.1
dt    = 0.1

xbig = linspace(0., 1., N+1)
x    = xbig[1:N]
dx   = 1./N

beta = alpha*dt/(2.*dx**2)

t0  = 0.
t1  = 1.5
phi0 = sin(pi*x)
phi  = phi0

A, f = makef(phi, beta)
for t in r_[t0:t1 + dt:dt]:
    if (t == 0.0 or t == .5 or t == 1.0 or t == 1.5):
        xexact = linspace(0., 1., 1000)
        uexact = multiply(sin(pi*xexact), exp(-alpha*pi**2*t))
        phibig = hstack((0., phi, 0.))
        annotate("t=%.2f" % t, xy = (x[x.shape[0]//2], phi.max()), xytext = (0., 5.), textcoords = "offset points")
        if t == 0.:
            plot(xbig, phibig, "ko", mfc = "none", label = "N=%d" % N)
            plot(xexact, uexact, "k", label = "exact")
        else:
            plot(xexact, uexact, "k")
            plot(xbig, phibig, "ko", mfc = "none")

    A, f = makef(phi, beta)
    phi  = tri_diag(A, f)

legend(loc = 0)
grid("on")
show()
