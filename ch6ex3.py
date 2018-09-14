# Example 6.3
# Calculation of the discrete Sine and Cosine transform
from pylab import*
from dtrans import*

# Define the function
f = lambda x: np.power((x/pi), 2)

# Sine and Cosine transform
N = 16
X = asmatrix(r_[0.:pi+pi/N:pi/N])
k = r_[0:N+1]
cos_coeff = dct1(f(X).T)
sin_coeff = dst1(f(X).T)


# Generates the plots
plot(k, abs(cos_coeff), 'k-^', lw = 2., label = "Cosine Transform")
plot(k, abs(sin_coeff), 'k:o', lw = 2., label = "Sine Transform")
grid("on")
xlim([-0.5, 16.5]); ylim([-0.1, 0.5])
legend(loc = 0)
show()
