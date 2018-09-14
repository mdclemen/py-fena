# Example 6.7
# Initial Boundary Value Problem
from pylab import*
from scipy.integrate import solve_ivp
# Initialization
clf()

N = 32
dX = 1./N

X = r_[0.:1.:dX]

# Initial Conditions

u_0 = zeros(N)
for j in range(len(X)):
    if X[j] < 0.4:
        u_0[j] = 1. - 25.*(X[j] - 0.2)**2

# Solve the ODE

k     = 2.*pi*r_[0.:N//2]
uk_0  = fft(u_0)/N
uk_0  = uk_0[0:N//2]

uk2 = hstack((uk_0, 0., conj(flipud(uk_0[1:N//2]))))
u = ifft(uk2*N)

fun = lambda t, y: -multiply(1j*k + 0.05*np.power(k, 2), y)
dt = .006
r = solve_ivp(fun, (0., .75), uk_0, method = "RK45", t_eval = r_[0:.75+dt:dt])

# Plot the solution (figure 6.9)

uk_sol = zeros([N//2, 1, 3], dtype=complex64)
u_sol  = zeros([N, 1, 3], dtype=complex64)
uk_sol[:,:,0] = r.y.T[41,:].reshape(N//2, 1)
uk_sol[:,:,1] = r.y.T[82,:].reshape(N//2, 1)
uk_sol[:,:,2] = r.y.T[124,:].reshape(N//2, 1)

for k in range(3):
    uk = hstack((uk_sol[:,:,k].flat, 0., conj(flipud(uk_sol[1:N//2,0,k]))))
    u_sol[:,:,k] = ifft(uk*N).reshape(uk.shape[0], 1)

figure(1)
plot(hstack((X, X[N-1]+dX)), hstack((u_0, u_0[0])), "k-o", lw = 2.,  label = "t=0")
plot(hstack((X, X[N-1]+dX)), hstack((u_sol[:,0,0], u_sol[0,0,0])), "k--o", lw = 2., label = "t=0.25")
plot(hstack((X, X[N-1]+dX)), hstack((u_sol[:,0,1], u_sol[0,0,1])), "k:o", lw = 2., label = "t=0.5")
plot(hstack((X, X[N-1]+dX)), hstack((u_sol[:,0,2], u_sol[0,0,2])), "k-.x", lw = 2., label = "t=0.75")
xlabel("x"); ylabel("u(x,t)")
 
legend(loc = 0); grid("on")
     
# Reproduce figure 6.8

figure(2)
h = r_[0:8e-3+5e-4:5e-4]
for n in [1., 5., 8., 11., 13., 14., 15.]:
    k = 2.*pi*n
    lam = -(1j*k + 0.05*k**2)
    lh = lam * h
    
    sigma = abs(1. + lh + (np.power(lh, 2))/2. + (np.power(lh, 3))/6. + (np.power(lh, 4))/24.)
    plot(h, sigma, "k", lw = 2.)
    xlabel("h"); ylabel("|sigma|")

show()
