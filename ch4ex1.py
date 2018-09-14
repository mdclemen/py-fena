# Chapter 4, Examples 1 & 2
# Explicit and Implicit Euler
from pylab import*

# Cleanup
clf()

# Analytic solution
f = lambda x: exp(-x/2.);
t_true = linspace(0.,20.,50);
y_true = f(t_true);

# Numerical solution for h=1
y1 = ones_like(y_true)
h = 1.
t1 = zeros_like(t_true)
for i in range(len(y1) - 1):
    y1[i+1] = y1[i] - .5*h*y1[i]
    t1[i+1] = t1[i] + h

# Numerical solution for h=4.2
y2 = ones_like(y_true)
h = 4.2
t2 = zeros_like(t_true)
for i in range(len(y2) - 1):
    y2[i+1] = y2[i] - .5*h*y2[i]
    t2[i+1] = t2[i] + h

# Plot explicit solutions
figure(1)
plot(t_true, y_true, label = "Exact")
plot(t1, y1, "k.--", ms = 10, label = "Explicit Euler, h=1")
plot(t2, y2, "ro--", ms = 10, mfc = "none", label = "Explicit Euler, h=4.2")
xlim([0., 20.]); ylim([-1.5, 2.])
xlabel("t")
ylabel("y(t)")
legend(loc = "upper left")

# Implicit solution for h=1
y3 = ones_like(y_true)
h = 1.
t3 = zeros_like(t_true)
for i in range(len(y3) - 1):
    y3[i+1] = y3[i]/(1.+h/2.)
    t3[i+1] = t3[i] + h

# Implicit solution for h=4.2
y4 = ones_like(y_true)
h = 4.2;
t4 = zeros_like(t_true)
for i in range(len(y4) - 1):
    y4[i+1] = y4[i]/(1.+h/2.)
    t4[i+1] = t4[i] + h

# Plot implicit solutions
figure(2)
plot(t_true, y_true, label = "Exact")
plot(t3, y3, "k.--", ms = 10, label = "Implicit Euler, h=1")
plot(t4, y4, "rs--", ms = 10, mfc = "none", label = "Implicit Euler, h=4.2")
xlim([0., 20.]); ylim([-.1, 1.1])
xlabel("t")
ylabel("y(t)")
legend(loc = "upper right")
show()
