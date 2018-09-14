# Example 6.1
from pylab import*
# Initialization
clf()

N = 16
X = r_[0. : 32.*pi/16. : 2*pi/16.]
Y = zeros(size(X))

# Define the function (step down)
for j in range(size(X)):
    if X[j] < (pi-0.001):
        Y[j] = 1.
    else:
        Y[j] = -1.

# Coefficients of the FT
YY = fft(Y)/N
Z = YY[0:N//2 + 1]
print("Z =\n", Z)

# Generate the plot

x = r_[-8:9:1]
Z = hstack((flipud(Z[1::]), Z))

plot(x, abs(Z), "k.")
grid("on")
show()
