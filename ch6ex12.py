# Example 6.12
# Calculation of Derivatives Using Chebyshev Derivative Matrix Operator
from pylab import*
# Initialization
clf()

N = 5

u = lambda x: 4.*multiply(np.power(x, 2), multiply(1. - np.power(x, 2), exp(-x/2.)))
j = r_[0:N+1]
X = cos(j*pi/N)
U = u(X)

D = asmatrix(zeros([N+1,N+1]))

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

# Display the results

print("X =\n", X)
print("U =\n", U)
print("D =\n", D)
print("D*U =\n", dot(D, U))
