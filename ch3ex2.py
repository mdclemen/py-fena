# Example 3.2
#  Romberg integration
from pylab import*

# I_1, I_2, I_3, I_4
I = zeros([4,1])
h = zeros([4,1])

for i in range(4):
    N    = 2**(i+1)
    h[i] = (pi - 1.)/N
    X    = r_[1.:pi + h[i]:h[i]]
    f    = sin(X)/(2.*np.power(X, 3))
    
    for j in range(1,N):
        I[i,:] = I[i,:] + f[j]
    I[i] = I[i] + 1./2. * (f[N] + f[0])
    I[i] = I[i]*h[i]

# I_12, I_23, I_24

I_12 = (4.*I[1] - I[0])/3.
I_23 = (4.*I[2] - I[1])/3.
I_34 = (4.*I[3] - I[2])/3.

# I_123, I_234

I_123 = (16.*I_23 - I_12)/15.
I_234 = (16.*I_34 - I_23)/15.

# I_1234

print("I_1234 = %.15f" % ((64.*I_234 - I_123)/63.))
