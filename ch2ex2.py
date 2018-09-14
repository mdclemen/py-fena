# Example 2.2
# Pade Differentiation Using a Lower Order Boundry Scheme
from pylab import*
# Setup
clf()

# function and derivative
f = lambda x: sin(5.*x)
df = lambda x: 5.*cos(5.*x)
# grid
N = 15
H = linspace(0. ,3., N)
h = 3./(N-1)

# sample function and derivative
func = f(H)
deriv = df(H)

# Construct difference scheme
# Here we are using a fourth-order Pade scheme with third order boundary
# schemes given by (2.17) in the text.

A = zeros([N,N])
b = zeros(N)
A[0,0:2] = array([1., 2.])
A[N-1,N-2:N] = array([2., 1.])
b[0] = -5.*func[0]/2. + 2.*func[1] + func[2]/2.
b[N-1] = 5.*func[N-1]/2. - 2.*func[N-2] - func[N-3]/2.
for i in range(1, N-1):
    A[i,i-1:i+2] = array([1., 4., 1.])
    b[i] = 3.*(func[i+1] - func[i-1])
b = b/h

# Compute approximate derivative
approx_deriv = linalg.solve(A, b)

# Generate Plots
figure(1); clf()
plot(H, deriv, "k.-", ms = 10, label = "Exact Derivative")
plot(H, approx_deriv, "--", label = "Computed Derivative")
xlabel("x", fontsize = 14)
ylabel("Derivative", fontsize = 14)
title("4th Order Pade Differentiation", fontsize = 14)
legend(loc = 0)
grid("on")
show()
