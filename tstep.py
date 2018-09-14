from pylab import*
def explicit_euler(f, tspan, y0, N):
# Explicit Euler

    m = size(y0)
    t = linspace(tspan[0], tspan[1], N+1)
    y = asmatrix(zeros([m,N+1], dtype = "complex"))
    h = (tspan[1] - tspan[0])/N
    y[:,0] = y0

    for i in range(N):
        y[:,i+1] = y[:,i] + h*f(t[i], y[:,i])

    return t, y

def rk2(f, tspan, y0, N):
# Second-order Runge-Kutta

    m = size(y0)
    t = linspace(tspan[0], tspan[1], N+1)
    y = asmatrix(zeros([m,N+1], dtype = "complex"))
    h = (tspan[1] - tspan[0])/N
    y[:,0] = y0

    alpha = 1./2.

    for i in range(N):
        k1 = h*f(t[i], y[:,i])
        k2 = h*f(t[i] + alpha*h, y[:,i] + alpha*k1)
        y[:,i+1] = y[:,i] + (1. - 1./2./alpha)*k1 + k2/2./alpha

    return t, y

def rk4(f, tspan, y0, N):
# Fourth-order Runge-Kutta

    m = size(y0)
    t = linspace(tspan[0], tspan[1], N+1)
    y = asmatrix(zeros([m,N+1]), dtype = "complex")
    h = (tspan[1] - tspan[0])/N
    y[:,0] = y0

    for i in range(N):
        k1 = h*f(t[i], y[:,i])
        k2 = h*f(t[i] + h/2., y[:,i] + k1/2.)
        k3 = h*f(t[i] + h/2., y[:,i] + k2/2.)
        k4 = h*f(t[i] + h, y[:,i] + k3)
        y[:,i+1] = y[:,i] + k1/6. + (k2+k3)/3. + k4/6.

    return t, y
