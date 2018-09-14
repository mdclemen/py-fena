# Example 4.5
# Amplification Factor
from pylab import*
from tstep import*
# Setup
clf()

g = lambda t, y: 1j*y

# Generate Solutions
t, yee = explicit_euler(g, [0., 20.], 1., 100)
t, yrk = rk2(g, [0., 20.], 1., 100)
ytrue = exp(1j*t)

# Amplification
amp_true = abs(ytrue).max()
amp_ee   = abs(yee).max()
amp_rk   = abs(yrk).max()
print("amp_true = %.4f" % amp_true)
print("amp_ee = %.4f" % amp_ee)
print("amp_rk = %.4f" % amp_rk)

# Plot
figure(1)
plot(t, real(yee.flat), label = "Euler")
plot(t, real(yrk.flat), "r", label = "RK2")
plot(t, real(ytrue), "g", label = "True")
grid("on")
xlabel("t")
ylabel("y")
legend(loc = "upper left")
show()
