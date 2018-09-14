# Example 4.9
# Shooting to Solve the Blasius Boundary Layer:
# f''' + f*f'' = 0
from pylab import*
from tstep import rk4
# Clean and Setup
clf()

N = 1000
niter = 8

# Function
f = lambda t, y: asmatrix([[-y[0,0]*y[2,0]], [y[0,0]], [y[1,0]]])

# Initial Solves
t, y1 = rk4(f, [0., 10.], asmatrix([[1.], [0.], [0.]]), N)
t, y2 = rk4(f, [0., 10.], asmatrix([[.5], [0.], [0.]]), N)
f1a = .5
f1b = 1.
f2a = y2[1,-1]
f2b = y1[1,-1]

# Shooting
tol = 1e-18
for i in range(niter):
    f1new = f1a + (f1a - f1b)/(f2a - f2b)*(1. - f2a)
    t, y = rk4(f, [0., 10.], asmatrix([[f1new], [0.], [0.]]), N)
    print("iter:%2d, f1(0)=%17.15f, f2(10)=%17.15f" % (i, f1new ,y[1,-1]))
    if (abs(y[1,-1] - 1.) <= tol):
        break
    f1b = f1a
    f1a = f1new
    f2b = f2a
    f2a = y[1,-1]

# Generate Plots
plot(t, y[0,:].flat, label = "f_1")
plot(t, y[1,:].flat, label = "f_2")
plot(t, y[2,:].flat, label = "f_3")
grid("on")
xlim([0., 6.]); ylim([0., 2.])
legend(loc = 0)
xlabel(r"$\eta$", fontsize = 14)
show()
