# Example 4.4
# Linearization
from pylab import*
# Clean up
clf()

# Obtain true solution using sympy if available
try:
    from sympy import dsolve, Derivative, symbols
    y, t = symbols("y t", function = True)
    y = dsolve(Derivative(y(t), t) + y(t)*(1 - y(t)), 0, ics={y(0): .5}).rhs
    print("y(t) =", y)
    true_solution = y.evalf(subs = {t:1.})
except ImportError:
    # suppy the true soln
    true_solution = 1./(exp(1.) + 1.)

# Setup schemes
# Linearized Trapeziodal method
lin_trap = lambda y, h: y + h*y*(y - 1.)/(1. - h*(y - 1./2.))
full_trap = lambda y, h: (2./h + 1. - sqrt((2./h + 1.)**2 - 4.*(2./h*y + y*(y - 1.))))/2.

# Collect Data
H = logspace(0., 2.7, 10).round()
STEPS = zeros(len(H))
LIN_ERR = zeros(len(H))
FULL_ERR = zeros(len(H))

for i in range(len(H)):
    N = H[i]
    h = 1./N
    yl = 1./2.
    yf = 1./2.
    for j in range(int(N)):
        yl = lin_trap(yl, h)
        yf = full_trap(yf, h)
    STEPS[i] = N
    LIN_ERR[i] = abs(yl - true_solution)
    FULL_ERR[i] = abs(yf - true_solution)

# Plot
figure(1); clf()
loglog(STEPS, LIN_ERR, "o-", mfc = "none", label = "Trapeziodal")
loglog(STEPS, FULL_ERR, "s-", mfc = "none", label = "Linearized Trapeziodal")
grid("on")
xlabel("N -- Number of Steps", fontsize = 14)
ylabel("Error", fontsize = 14)
legend(loc = 0)
show()
