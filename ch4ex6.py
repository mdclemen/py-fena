# Example 4.6
# Runge-Kutta
from pylab import*
from tstep import*
# Setup
clf()

h = 0.15
w = 4.
Tmax = 6.

# Analyic solution
g = lambda x: cos(w*x)
T = linspace(0., 6., 100)
y_true = g(T)

# Function
f = lambda t, y: asmatrix([[y[1,0]], [-16.*y[0,0]]])

# Use Runge-Kutta code
t, y2 = rk2(f, [0., 6.], asmatrix([[1.], [0.]]), 40)
t, y4 = rk4(f, [0., 6.], asmatrix([[1.], [0.]]), 40)
yex = exp(4.*1j*t)

# Plot
figure(1)
plot(t , y2[0,:].flat, "--", label = "2nd O Runge-Kutta")
plot(t, y4[0,:].flat, "^:", mfc = "none", label = "4th O Runge-Kutta")
plot(T, y_true, label = "Exact")
grid("on")
xlim([0., 6.]); ylim([-2., 3.])
xlabel("t")
ylabel("y(t)")
legend(loc = "upper left")
show()
