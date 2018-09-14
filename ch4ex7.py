# Example 4.7
# Multi-Step Methods
# Same problem as example 4.3
from pylab import*
# Setup
clf()

h = 0.1
w = 4.
Tmax = 6.
N = int(Tmax / h)
t = linspace(0., 6., N+1)

# Analyic solution
g = lambda x: cos(w*x)
T = linspace(0., 6. , 100)
y_true = g(T)

# Function
f = lambda t, y: asmatrix([[y[1,0]], [-16.*y[0,0]]])

# Leapfrog
# initial values
Ylf = asmatrix(zeros([2,N+1]))
Ylf[:,0] = asmatrix([[1], [0]])
# take first step
Ylf[:,1] = Ylf[:,0] + h*f(t[0], Ylf[:,0])
# take remaining leapfrog steps
for i in range(1,N):
    Ylf[:,i+1] = Ylf[:,i-1] + 2.*h*f(t[i],Ylf[:,i])

# Adams-Bashforth
# initial values
Yab = asmatrix(zeros([2,N+1]))
Yab[:,0] = asmatrix([[1], [0]])
# take first step
Yab[:,1] = Yab[:,0] + h*f(t[0], Yab[:,0])
# take remaining leapfrog steps
for i in range(1,N):
    Yab[:,i+1] = Yab[:,i] + 3./2.*h*f(t[i],Yab[:,i]) - h/2.*f(t[i-1], Yab[:,i-1])

# Plot Generation
figure(1)
plot(t, Ylf[0,:].flat, "o:", mfc = "none", label = "Leapfrog")
plot(t, Yab[0,:].flat, "--", label = "Adams-Bashforth")
plot(T, y_true, label = "Exact")
xlim([0., 6.]); ylim([-2., 3.])
grid("on")
legend(loc = "upper left")
xlabel("t")
ylabel("y(t)")
show()
