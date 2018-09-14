# Example 4.8
# A Stiff System (Byrne and Hindmarsh)
# Here we use ode23s
from pylab import*
from scipy.integrate import solve_ivp
# Settings and Function
clf()

a = 1.5e-18
b = 2.5e-6
g = 2.1e-6
r = .6
s = .18
u = .016

# function
f = lambda t, y: array([-y[0]*(a*y[1] + b) + g, y[1]*(r*y[0] - s) + u*(1. + y[0])])

# Jacobian
J = lambda t, y: asmatrix([[-(a*y[1] + b), -a*y[0]], [r*y[1] + u, r*y[0] - s]])

# ODE Solver
# Use a scipy.integrate.solve_ivp with RK23 method (close to ode23s, stiff solver)
rode = solve_ivp(f, (0., 7.e5), array([-1., 0.]), method = "RK23", jac = J, rtol = 1.e-4, atol = 1.e-4*ones(2))

# Population Plot
figure(1)
plot(rode.t, rode.y[0,:].flat)
xlabel("Time")
ylabel("Population Inversion")
grid("on")

# Photon Density Plot
figure(2)
semilogy(rode.t, rode.y[1,:].flat)
xlabel("Time")
ylabel("Photon Density")
grid("on")
show()
