# Example 5.1
# Inhomogeneous heat equation
from pylab import*
from scipy.sparse import diags
# Setup
clf()
N = 21
L = 1.
alpha = 1
dx = L/(N-1)
X = asmatrix(linspace(0.,L,N)).T

# initial condition
T = sin(pi*X)

# Spatial derivative operator
A = asmatrix(diags(array([ones(N-3), -2.*ones(N-2), ones(N-3)]), [-1, 0, 1]).toarray())

# inhomogeneous term
f = lambda t: (pi**2-1)*exp(-t)*sin(pi*X[1:N-1,:])

# time advancement
g = lambda t,dt,T: T[1:N-1,:] + alpha*dt/(dx**2)*A*T[1:N-1,:] + dt*f(t)

# Stable run
dt = 0.001           # time step
t_final = 2.0        # final time
time = r_[0.:t_final+dt:dt]  # time array
pt = array([0.0, 0.5, 1.0, 1.5, 2.0])   # desired plot times
pn = len(pt)         # number of desired plots
pc = 0               # plot counter
rt = asmatrix(zeros([1, pn]))
S = asmatrix(zeros([N, pn]))   # solution storage

for t in time:
    # plot storage
    if (t >= pt[pc]):
        S[:,pc] = T
        rt[0,pc] = t
        pc += 1
        if (pc > pn):
            break

    # time advancement
    T[1:N-1,:] = g(t,dt,T)

# Plot stable run
figure(1)
linstyl = ["-", "--", ":", "-.", ".-"]
for i in range(pn):
    plot(X, S[:,i], "k%s" % linstyl[i], lw = 1., label = "t = %.1f" % rt[0,i])
xlabel("x", fontsize = 14)
ylabel("T(x)", fontsize = 14)
legend(loc = 0)
xlim([0., 1.]); ylim([0., 1.1])
ax = gca()
ax.set_xticks([0., .25, .5, .75, 1.])
ax.set_yticks([0., .25, .5, .75, 1.])
grid(axis="both")

# Unstable run
T = sin(pi*X)            # reset initial condition
dt = 0.0015              # time step
t_final = .3             # final time
time = r_[0.:t_final+dt:dt] # time array
pt = array([0.0, 0.15, 0.153, 0.158, 0.166]) # desired plot times
pn = len(pt)             # number of desired plots
pc = 0                   # plot counter
rt = asmatrix(zeros([1 ,pn]))
S = asmatrix(zeros([N, pn]))       #solution storage

for t in time:
    # plot storage
    try:
        if (t >= pt[pc]):
            S[:,pc] = T
            rt[0,pc] = t
            pc += 1
            if (pc > pn):
                break

        # time advancement
        T[1:N-1,:] = g(t,dt,T)
    except:
        pass
# Plot unstable run
figure(2)
for i in range(pn):
    plot(X, S[:,i], "k%s" % linstyl[i], lw = 1., label = "t = %.1f" % rt[0,i])
xlabel("x", fontsize = 14)
ylabel("T(x)", fontsize = 14)
legend(loc = 0)
xlim([0., 1.]); ylim([0., 1.3])
ax = gca()
ax.set_xticks([0., .25, .5, .75, 1.])
ax.set_yticks([0., .25, .5, .75, 1., 1.25])
grid(axis="both")
show()
