# Example 6.13
# Convection Equation with non constant coefficients
from pylab import*
# Initialization

clf()

N  = 16
h  = 0.001

u  = lambda x,t: sin(2.*pi*x*exp(-2.*t))
x  = r_[-1:1. + .01:.01]

j  = r_[0:N+1]
X  = cos(j*pi/N)
U  = u(X,0)

D  = zeros([N+1,N+1])

# Creation of the matrix operator

for j in range(N+1):
    for k in range(N+1):
        if j == k:
            if j == 0:
                D[j,k] = (2.*N**2 + 1.)/6.
            elif j == N:
                D[j,k] = -(2.*N**2 + 1.)/6.
            else:
                D[j,k] = -X[j]/(2.*(1. - X[j]**2))
        else:
            if (j == 0 or j == N):
                cj = 2
            else:
                cj = 1.
            if (k == 0 or k == N):
                ck = 2
            else:
                ck = 1.
            D[j,k] = cj*(-1.)**((j + 1) + (k + 1))/(ck*(X[j]-X[k]))

# Initial condition

plot(x, u(x,0), "k")
xlim([-1., 1.]); ylim([-1.5, 1.5])

# Iterate and plot

A  = eye(N+1) - 2.*h*dot(diag(X), D)

for t in r_[h:.6 + h:h]:
    U = dot(A, U)
  
    if t == 0.3:
        plot(x, u(x,t), "k--", X, U, "ko", mfc = "none")
    if t == 0.6:
        plot(x, u(x,t), "k--", X, U, "k*", mfc = "none")

grid("on")
show()    
