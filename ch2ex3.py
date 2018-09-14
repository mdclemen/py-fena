# Example 2.3
# Calculation of Derivatives on a Non-uniform Mesh
from pylab import*
# Setup
N = 18
J = r_[0:N+1]
y = .9*(2*J/N - 1.)
x = arctanh(y)

# Plot
figure(1); clf()
plot(y, x ,"*-", mfc = "none")
ylim([x.min(), y.max()]); ylim([x.min(), x.max()])
xlabel(r"$\xi$", fontsize = 14)
ylabel("x", fontsize = 14)
grid("on")

# Numbers
print(hstack((asmatrix(x).T, asmatrix(y).T)))
show()
