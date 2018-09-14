# Example 6.8
# Burgers Equations
from pylab import*
# Initialization
clf()

N  = 32
f  = lambda x: 10.*sin(x)
dX = 2.*pi/N
X  = r_[0.:2.*pi:dX]
dt = 0.005
t  = 0.

x  = r_[0.:2.*pi + .01:.01]
plot(x, f(x), "k", lw = 2., label = "t=0")

# Create Differentiation matrix :

D = zeros([N,N])
for l in range(N):
    for j in range(N):
        if l != j:
            D[l,j] = 0.5*(-1.)**(l-j)/tan(pi*float(l-j)/N)

# Iterate and plot :

u   = f(X)

for t in r_[dt:1. + dt:dt]:
    u = u + dt*(dot(dot(D, D), u) - dot(dot(diag(u), D), u))
    
    if t == .1:
        plot(X, u, "ko", mfc = "none", label = "t=0.1")
    elif t == .2:
        plot(X, u, "kx", mfc = "none",label = "t=0.2")
    elif t == .4:
        print("here")
        plot(X, u, "k*", mfc = "none", label = "t=0.4")
    elif t == .6:
        plot(X, u, "ks", mfc = "none", label = "t=0.6")
    
legend(loc = 0); grid("on")
show()
