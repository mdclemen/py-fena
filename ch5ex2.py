# Example 5.2
# Convection Equation
from pylab import*
from scipy.sparse import diags
# Setup
clf()
L = 1.
dx = 0.01
dt = 0.01
c = 1.

X = r_[0.:L+dx:dx]
N = len(X)
B = asmatrix(diags(array([-1.*ones(N-2), zeros(N-1), ones(N-2)]), [-1, 0, 1]).toarray())

B[-1,-1] = 2.
B[-1,-2] = -2.

u = asmatrix(exp(-200.*(X-0.25)**2)).T
u[0] = 0.

# System of ODEs
f = lambda u: vstack([0., -c/2./dx*B*u[1::]])

# Explicit Euler
t_final = .4
time = r_[0.:t_final+dt:dt]
pt = array([0., .12, .25, .375])
pn = len(pt)
pc = 0
rt = asmatrix(zeros([1, pn]))
S = asmatrix(zeros([N, pn]))

for t in time:
    # plot storage
    try:
        if (t >= pt[pc]):
            S[:,pc] = u
            rt[0,pc] = t
            pc += 1
            if (pc > pn):
                break
    except:
        pass
        # time advancement
    u = u + dt*f(u)

# Plot unstable figure
figure(1)
linstyl = ["-", "--", ":", "-.", ".-"]
for i in range(pn):
    plot(X, S[:,i], "k%s" % linstyl[i], lw = 1., label = "t = %.2f" % rt[0,i])
xlim([0., .75]); ylim([-4., 4.])
ax = gca()
ax.set_xticks([0., .25, .5, .75])
ax.set_yticks([-4., -2., 0., 2., 4.])
grid(axis="both")
xlabel("x", fontsize = 14)
ylabel("u(x)", fontsize = 14)
legend(loc = "upper left")

# Fourth-order Runge-Kutta

# reset initial condition
u = asmatrix(exp(-200.*(X-0.25)**2)).T
u[0] = 0.

t_final = .8
time = r_[0:t_final+dt:dt]
pt = array([0., .25, .50, .75])
pn = len(pt)
pc = 0
rt = asmatrix(zeros([1, pn]))
S = asmatrix(zeros([N, pn]))

for t in time:
    # plot storage
    try:
        if (t >= pt[pc]):
            S[:,pc] = u
            rt[0,pc] = t
            pc += 1
            if (pc > pn):
                break
    except:
        pass
    # time advancement
    u1 = dt*f(u)
    u2 = dt*f(u+u1/2.)
    u3 = dt*f(u+u2/2.)
    u4 = dt*f(u+u3)
    u = u + (u1 + 2.*u2 + 2.*u3 + u4)/6.

# Plot stable figure
figure(2)
for i in range(pn):
    plot(X, S[:,i], "k%s" % linstyl[i], lw = 1., label = "t = %.2f" % rt[0,i])
xlim([0., 1.]); ylim([-.1, 1.5])
ax = gca()
ax.set_xticks([0., .25, .5, .75, 1.])
ax.set_yticks([0., .25, .5, .75, 1., 1.25, 1.5])
grid("on")
xlabel("x", fontsize = 14)
ylabel('u(x)', fontsize = 14)
legend(loc = "upper left")
show()
