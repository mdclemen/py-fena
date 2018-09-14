# Example 3.1
# Numerical integration of int(sin(x)/(2*x^3),0,pi)
from pylab import*
# Initialization
clf()

N = 8

h = (pi - 1.)/N  # uniform spacing

X = r_[1.:pi + h:h]
f = sin(X)/(2.*np.power(X, 3))

# Resolution

# A : Trapezoidal method
# B : Simpson
# C : End-corrected Trapezoidal
#
# err_A, err_B, err_C : % of error with 'exact' sol

I = 0.1985572988
A = 0         # trapezoidal
B = 0         # quadrature
C = 0         # correct end

#  Trapezoidal Rule
#--------------------
for i in range(1,N):
    A = A + f[i]
A = A + 1./2. * (f[N] + f[0])
A = A*h
err_A = abs(I - A)/I * 100.

#  Simpson (quad)
#--------------------
for i in range(N//2-1):
    B = B + 4.*f[2*i+1] + 2.*f[2*i+2]
B = B + f[0] + f[N] + 4.*f[N-1]
B = B*h/3.
err_B = abs(I - B)/I * 100.

#  Trap. End Correct
#--------------------

df_a = (cos(1.)-3.*sin(1.))/2.
df_b = -1./(2.*pi**3)

C = A - (h**2)/12. * (df_b - df_a)
err_C = abs(I - C)/I * 100.

# Display results

print("N = %d         Result     Error" % N)
print("Trapezoidal    %.5f       %.5f" % (A, err_A))
print("Simpson        %.5f       %.5f" % (B, err_B))
print("Trap end corr  %.5f       %.5f" % (C, err_C))
