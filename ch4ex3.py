# Example 4.3
# A second order equation
from pylab import*
# Setup
clf()

h = 0.15
w = 4.
Tmax = 6.

# Analyic solution
f = lambda x: cos(w*x)
T = linspace(0.,6.,100)
y_true = f(T)

# Explicit Euler
Te = r_[0.:6.+h:h]
Ye = asmatrix(zeros([2, size(Te)]))
Ye[:,0] = asmatrix([[1.],[0]])
A = asmatrix([[1., h], [-w**2*h, 1.]])
for i in range(size(Te) - 1):
    Ye[:,i+1] = A*Ye[:,i]

# Implicit Euler
Yi = asmatrix(zeros([2, size(Te)]))
Yi[:,0] = asmatrix([[1.],[0]])
B = asmatrix([[1., -h], [w**2*h, 1.]])
for i in range(size(Te) - 1):
    Yi[:,i+1] = linalg.solve(B, Yi[:,i])

# Trapezoidal
Yt = asmatrix(zeros([2, size(Te)]))
Yt[:,0] = asmatrix([[1.],[0]])
C = asmatrix([[1., -h/2.], [w**2*h/2., 1.]])
D = asmatrix([[1., h/2.], [-w**2*h/2., 1.]])
for i in range(size(Te) - 1):
    Yt[:,i+1] = linalg.solve(C, D)*Yt[:,i]

# Plot generation
figure(1)
plot(T, y_true, "k", label = "Exact")
plot(Te,Ye[0,:].flat,"r:", label = "Explicit Euler")
plot(Te,Yi[0,:].flat,"g-.", label ="Implicit Euler")
plot(Te,Yt[0,:].flat,"--", label = "Trapezoidal")
xlim([0., 6.]); ylim([-2.6, 2.6])
xlabel("t")
ylabel("y(t)")
legend(loc = 0)
show()
