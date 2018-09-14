# Example 2.1
# Accuracy of Finite Difference Schemes
from pylab import*
# Define functions

# Function
f = lambda x: sin(x)/x**3
# Analytic first derivative
df = lambda x: cos(x)/x**3 - 3.*sin(x)/x**4
# Forward difference
forward = lambda x, h: (f(x + h) - f(x))/h
# Central difference
central = lambda x, h: ( f(x + h) - f(x - h))/2./h
# Fourth order difference
fourth = lambda x, h: (f(x - 2.*h) - 8.*f(x - h) + 8.*f(x + h) - f(x + 2.*h) )/12./h

# Evaluate error in finite differences
x = 4

H_for = logspace(-5., 0., 20)
H_cen = logspace(-5., 0., 10)
H_fth = logspace(-2.5, 0., 10)
for_err = zeros(size(H_for))
cen_err = zeros(size(H_cen))
fth_err = zeros(size(H_fth))

for i in range(size(H_for)):
    for_err[i] = abs(forward(x, H_for[i]) - df(x))

for i in range(size(H_cen)):
    cen_err[i] = abs(central(x, H_cen[i]) - df(x))

for i in range(size(H_fth)):
    fth_err[i] = abs(fourth(x, H_fth[i]) - df(x))

# Generate plot
figure(1); clf()
loglog(H_for, for_err, "o-", mfc = "none", label = "1st order")
loglog(H_cen, cen_err, "ko-", mfc = "none", label = "2nd order")
loglog(H_fth, fth_err, "ro-", mfc = "none", label = "3rd order")
xlabel("h -- grid spacing", fontsize = 12)
ylabel("error", fontsize = 12)
legend(loc = "upper left")
grid("on")

# Compute slopes
print("first_order = %.6f" % mean(diff(log(for_err))/diff(log(H_for))))
print("second_order = %.6f" % mean(diff(log(cen_err))/diff(log(H_cen))))
print("fourth_order = %.6f" %  mean(diff(log(fth_err))/diff(log(H_fth))))
show()
